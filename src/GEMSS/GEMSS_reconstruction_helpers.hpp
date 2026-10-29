
/**
 * @file GEMSS_reconstruction_helpers.hpp
 * @brief Helper functions and the iterative solver of the GEMSS reconstruction.
 *
 * Provides peak detection and filtering, sphere table management, the summed field and MultisphereSolver.
 *
 * @author Arash Moradian
 * @date 2026-03-09
 */

#ifndef GEMSS_RECONSTRUCTION_HELPERS_HPP
#define GEMSS_RECONSTRUCTION_HELPERS_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <Eigen/Core>
#include <stdexcept>
#include <string>
#ifdef HAVE_OPENMP
    #include <omp.h>
#endif
#include "GEMSS_datatypes.hpp"
#include "GEMSS_config.hpp"
#include "GEMSS_io.hpp"
#include "GEMSS_voxel_processing.hpp"
#include "GEMSS_edt_incremental.hpp"

namespace GEMSS {

namespace peak_detail {

/// Search window of a peak test: the offsets of the (2w+1)^3 neighbourhood, nearest first.
struct Window {
    int w = 0;                             ///< half width (search_window)
    std::vector<Eigen::Vector3i> offsets;  ///< (dx, dy, dz), sorted by squared distance
    std::vector<long> linear;              ///< the same offsets as linear index steps
};

/**
 * @brief Builds the search window, ordered nearest-first (stable sort, so ties keep a fixed order).
 * @param search_window Half width w of the window.
 * @param ny, nz Grid size along y and z (for the linear index steps; nz == 1 is a 2D grid without z-neighbours).
 */
inline Window make_window(int search_window, int ny, int nz) {
    Window win;
    win.w = std::max(search_window, 0);
    const int wz = (nz == 1) ? 0 : win.w;
    const long sy = static_cast<long>(ny) * nz, sz = nz;
    for (int dx = -win.w; dx <= win.w; ++dx)
        for (int dy = -win.w; dy <= win.w; ++dy)
            for (int dz = -wz; dz <= wz; ++dz)
                if (dx || dy || dz) win.offsets.emplace_back(dx, dy, dz);
    std::stable_sort(win.offsets.begin(), win.offsets.end(), [](const Eigen::Vector3i& a, const Eigen::Vector3i& b) {
        return a.squaredNorm() < b.squaredNorm();
    });
    for (const Eigen::Vector3i& o : win.offsets) win.linear.push_back(o.x() * sy + o.y() * sz + o.z());
    return win;
}

/**
 * @brief Peak test of one voxel.
 *
 * A voxel is a peak if its distance is at least min_radius_vox, its field value exceeds 1, it is outside the
 * overlap mask, and no in-bounds voxel of its window has a larger field value. On a smooth field a non-peak has
 * a larger immediate neighbour, so the nearest-first window rejects it after about two reads.
 *
 * @param F Field accessor, F[i] with i the linear index (z fastest).
 * @param D Original distance field (radius values).
 * @param M Overlap mask.
 * @param win Search window from make_window().
 * @param x, y, z Voxel coordinates; i Linear index of the voxel.
 * @param nx, ny, nz Grid size.
 * @param config Uses min_radius_vox.
 * @param radius Set to D[i] (also when the test fails early).
 * @return True if the voxel is a peak.
 */
template <typename Field>
inline bool is_peak(const Field& F, const float* D, const uint8_t* M, const Window& win,
                    int x, int y, int z, long i, int nx, int ny, int nz,
                    const MultisphereConfig& config, float& radius) {
    radius = D[i];
    if (radius < config.min_radius_vox) return false;
    const float value = F[i];
    if (value <= 1.0f) return false;
    if (M[i] > 0) return false;
    if (radius > win.w) {
        // D > w: the whole window lies inside the grid (outside the grid is background), so no bound checks
        for (const long step : win.linear) {
            if (F[i + step] > value) return false;
        }
        return true;
    }
    // Only reachable when min_radius_vox <= search_window: skip neighbours outside the grid
    const long sy = static_cast<long>(ny) * nz, sz = nz;
    for (const Eigen::Vector3i& o : win.offsets) {
        const int ix = x + o.x(), iy = y + o.y(), iz = z + o.z();
        if (ix < 0 || ix >= nx || iy < 0 || iy >= ny || iz < 0 || iz >= nz) continue;
        if (F[ix * sy + iy * sz + iz] > value) return false;
    }
    return true;
}

/// Peak order: radius descending, then x, y, z ascending.
inline bool peak_order(const Eigen::Vector4f& a, const Eigen::Vector4f& b) {
    const float epsilon = 1e-5f;
    if (std::abs(a[3] - b[3]) > epsilon) return a[3] > b[3];
    if (a[0] != b[0]) return a[0] < b[0];
    if (a[1] != b[1]) return a[1] < b[1];
    return a[2] < b[2];
}

/// Sorts the peaks with peak_order() and returns them as an Nx4 matrix (x, y, z, radius).
inline Eigen::MatrixX4f sorted_matrix(std::vector<Eigen::Vector4f>& peaks) {
    std::sort(peaks.begin(), peaks.end(), peak_order);
    Eigen::MatrixX4f result(peaks.size(), 4);
    for (size_t i = 0; i < peaks.size(); ++i) result.row(i) = peaks[i];
    return result;
}

} // namespace peak_detail

/**
 * @brief Detects local maxima (peaks) of a field in a 3D volume.
 *
 * Runs peak_detail::is_peak() on every voxel. The result does not depend on the thread count.
 *
 * @param F Field accessor, F[i] with i the linear index (z fastest); same shape as original_distance.
 * @param original_distance Original distance field (radius values).
 * @param overlap_mask Voxels excluded from being peaks (near accepted spheres).
 * @param config Uses search_window, min_radius_vox and radius_offset_vox.
 * @return Nx4 matrix (x, y, z, radius) of the peaks, radius = distance + radius_offset_vox, sorted by
 *         radius (descending), then x, y, z.
 */
template <typename Field>
inline Eigen::MatrixX4f peak_local_max_3d(
    const Field& F,
    const VoxelGrid<float>& original_distance,
    const VoxelGrid<uint8_t>& overlap_mask,
    const MultisphereConfig& config)
{
    const int nx = static_cast<int>(original_distance.nx());
    const int ny = static_cast<int>(original_distance.ny());
    const int nz = static_cast<int>(original_distance.nz());
    const long sy = static_cast<long>(ny) * nz, sz = nz;
    const peak_detail::Window win = peak_detail::make_window(config.search_window, ny, nz);
    const float* D = original_distance.data.data();
    const uint8_t* M = overlap_mask.data.data();

    std::vector<Eigen::Vector4f> peaks;
    #pragma omp parallel
    {
        std::vector<Eigen::Vector4f> local_peaks;
        #pragma omp for collapse(2) schedule(dynamic, 16)
        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < ny; ++y) {
                const long row = x * sy + y * sz;
                for (int z = 0; z < nz; ++z) {
                    float radius;
                    if (peak_detail::is_peak(F, D, M, win, x, y, z, row + z, nx, ny, nz, config, radius)) {
                        local_peaks.emplace_back(static_cast<float>(x), static_cast<float>(y), static_cast<float>(z),
                                                 radius + config.radius_offset_vox);
                    }
                }
            }
        }
        #pragma omp critical
        peaks.insert(peaks.end(), local_peaks.begin(), local_peaks.end());
    }
    return peak_detail::sorted_matrix(peaks);
}

/**
 * @brief Incremental peak detection: re-tests only the voxels of the marked blocks and keeps the previous peaks
 *        outside them that are still outside the overlap mask.
 *
 * Gives the same result as peak_local_max_3d() if, since the previous search, the overlap mask only grew and the
 * marked blocks contain every voxel within the search window of a voxel whose field value changed. Outside the
 * marked blocks a voxel's value and its whole window are then unchanged, so its peak status can only be lost by
 * masking.
 *
 * @param F Field accessor, F[i] with i the linear index (z fastest); same shape as original_distance.
 * @param original_distance Original distance field (radius values).
 * @param overlap_mask Voxels excluded from being peaks.
 * @param config Uses search_window, min_radius_vox and radius_offset_vox.
 * @param test_blocks Flags of the blocks to re-test, indexed (bx * nby + by) * nbz + bz.
 * @param block Block edge length (voxels).
 * @param previous Peaks (x, y, z, radius) of the previous search; replaced by the new peaks.
 * @return The new peaks, in the same form and order as peak_local_max_3d().
 */
template <typename Field>
inline Eigen::MatrixX4f peak_local_max_3d_incremental(
    const Field& F,
    const VoxelGrid<float>& original_distance,
    const VoxelGrid<uint8_t>& overlap_mask,
    const MultisphereConfig& config,
    const std::vector<uint8_t>& test_blocks,
    int block,
    std::vector<Eigen::Vector4f>& previous)
{
    const int nx = static_cast<int>(original_distance.nx());
    const int ny = static_cast<int>(original_distance.ny());
    const int nz = static_cast<int>(original_distance.nz());
    const long sy = static_cast<long>(ny) * nz, sz = nz;
    const int nby = (ny + block - 1) / block, nbz = (nz + block - 1) / block;
    const peak_detail::Window win = peak_detail::make_window(config.search_window, ny, nz);
    const float* D = original_distance.data.data();
    const uint8_t* M = overlap_mask.data.data();

    // Keep the previous peaks outside the marked blocks that are still unmasked.
    std::vector<Eigen::Vector4f> peaks;
    for (const Eigen::Vector4f& p : previous) {
        const int x = static_cast<int>(p[0]), y = static_cast<int>(p[1]), z = static_cast<int>(p[2]);
        const long b = (static_cast<long>(x / block) * nby + y / block) * nbz + z / block;
        if (!test_blocks[b] && M[x * sy + y * sz + z] == 0) peaks.push_back(p);
    }

    // Re-test every voxel of the marked blocks.
    std::vector<int> blocks;
    for (size_t b = 0; b < test_blocks.size(); ++b) {
        if (test_blocks[b]) blocks.push_back(static_cast<int>(b));
    }
    #pragma omp parallel
    {
        std::vector<Eigen::Vector4f> local_peaks;
        #pragma omp for schedule(dynamic, 1)
        for (size_t k = 0; k < blocks.size(); ++k) {
            const int bx = blocks[k] / (nby * nbz), by = (blocks[k] / nbz) % nby, bz = blocks[k] % nbz;
            const int x_end = std::min(nx, (bx + 1) * block);
            const int y_end = std::min(ny, (by + 1) * block);
            const int z_end = std::min(nz, (bz + 1) * block);
            for (int x = bx * block; x < x_end; ++x) {
                for (int y = by * block; y < y_end; ++y) {
                    const long row = x * sy + y * sz;
                    for (int z = bz * block; z < z_end; ++z) {
                        float radius;
                        if (peak_detail::is_peak(F, D, M, win, x, y, z, row + z, nx, ny, nz, config, radius)) {
                            local_peaks.emplace_back(static_cast<float>(x), static_cast<float>(y), static_cast<float>(z),
                                                     radius + config.radius_offset_vox);
                        }
                    }
                }
            }
        }
        #pragma omp critical
        peaks.insert(peaks.end(), local_peaks.begin(), local_peaks.end());
    }
    Eigen::MatrixX4f result = peak_detail::sorted_matrix(peaks);
    previous.swap(peaks);
    return result;
}

/**
 * @brief Detects local maxima of a stored field (see the templated peak_local_max_3d()).
 * @param field Field to search for peaks.
 * @param original_distance Original distance field (radius values).
 * @param overlap_mask Voxels excluded from being peaks.
 * @param config Uses search_window, min_radius_vox and radius_offset_vox.
 * @return Nx4 matrix (x, y, z, radius) of the peaks.
 */
inline Eigen::MatrixX4f peak_local_max_3d(
    const VoxelGrid<float>& field,
    const VoxelGrid<float>& original_distance,
    const VoxelGrid<uint8_t>& overlap_mask,
    const MultisphereConfig& config)
{
    return peak_local_max_3d(GridField{field.data.data()}, original_distance, overlap_mask, config);
}

/**
 * @brief Filters peaks based on minimum center distance and applies sub-voxel shifts.
 * @param peaks Candidate peaks (Nx4 matrix), in the order of peak_local_max_3d().
 * @param sphere_table History table of accepted spheres.
 * @param F Field accessor for the sub-voxel shift, F[i] with i the linear index (z fastest).
 * @param nx, ny, nz Grid size of F.
 * @param config Uses min_center_distance_rel.
 * @return Filtered peaks as Nx4 matrix (x, y, z, radius).
 */
template <typename Field>
inline Eigen::MatrixX4f filter_and_shift_peaks(
    const Eigen::MatrixX4f& peaks,
    const Eigen::MatrixX4f& sphere_table,
    const Field& F, int nx, int ny, int nz,
    const MultisphereConfig& config)
{
    if (peaks.rows() == 0) return peaks;
    const float min_dist_sq = (float)(config.min_center_distance_rel * config.min_center_distance_rel);
    std::vector<Eigen::Vector4f> accepted_spheres;
    accepted_spheres.reserve(peaks.rows());

    for (int i = 0; i < peaks.rows(); ++i) {
        float px = peaks(i, 0);
        float py = peaks(i, 1);
        float pz = peaks(i, 2);

        bool collision = false;

        // Check self-collision
        if (!collision) {
            for (const auto& sphere : accepted_spheres) {
                float dx = px - sphere[0];
                float dy = py - sphere[1];
                float dz = pz - sphere[2];
                if ((dx*dx + dy*dy + dz*dz) < min_dist_sq * sphere[3]) { // Scale by sqrt (radius) to allow closer placement of smaller spheres
                    collision = true;
                    break;
                }
            }
        }
        if (!collision) {
            Eigen::Vector3i coord_int = peaks.row(i).head<3>().cast<int>();
            Eigen::Vector3f direction = shift_voxel_center(F, nx, ny, nz, coord_int, 0.25f);
            Eigen::Vector3f shift_vec = 0.5f * direction;
            Eigen::Vector3f true_center = coord_int.cast<float>() + shift_vec;
            float true_radius = peaks(i, 3) + shift_vec.norm();
            accepted_spheres.push_back(Eigen::Vector4f(
                true_center.x(),
                true_center.y(),
                true_center.z(),
                true_radius
            ));
        }
    }

    Eigen::MatrixX4f result(accepted_spheres.size(), 4);
    for (size_t k = 0; k < accepted_spheres.size(); ++k) {
        result.row(k) = accepted_spheres[k];
    }
    return result;
}

/**
 * @brief Filters peaks based on minimum center distance and applies sub-voxel shifts on a stored field.
 * @param peaks Candidate peaks (Nx4 matrix).
 * @param sphere_table History table of accepted spheres.
 * @param distance_field Field for the sub-voxel shift.
 * @param config Uses min_center_distance_rel.
 * @return Filtered peaks as Nx4 matrix (x, y, z, radius).
 */
inline Eigen::MatrixX4f filter_and_shift_peaks(
    const Eigen::MatrixX4f& peaks,
    const Eigen::MatrixX4f& sphere_table,
    VoxelGrid<float>& distance_field,
    MultisphereConfig config)
{
    return filter_and_shift_peaks(peaks, sphere_table, GridField{distance_field.data.data()},
                                  static_cast<int>(distance_field.nx()), static_cast<int>(distance_field.ny()),
                                  static_cast<int>(distance_field.nz()), config);
}

/**
 * @brief Filters the sphere table to keep only the largest connected network.
 * Uses an OpenMP-parallelized Dense Grid Spatial Hash with dynamic search extents.
 * @param sphere_table Input sphere table (Nx4 matrix).
 * @return Filtered sphere table containing only the primary connected component.
 */
inline Eigen::MatrixX4f filter_largest_sphere_network(const Eigen::MatrixX4f& sphere_table) {
    int n = sphere_table.rows();
    if (n <= 1) return sphere_table;

    // Disjoint-Set (Union-Find)
    struct UnionFind {
        std::vector<int> parent;
        std::vector<float> volume;

        UnionFind(int size) : parent(size), volume(size, 0.0f) {
            for (int i = 0; i < size; ++i) parent[i] = i;
        }
        int find(int i) {
            int root = i;
            while (root != parent[root]) root = parent[root];
            int curr = i;
            while (curr != root) {
                int nxt = parent[curr];
                parent[curr] = root;
                curr = nxt;
            }
            return root;
        }
        void unite(int i, int j) {
            int root_i = find(i);
            int root_j = find(j);
            if (root_i != root_j) {
                parent[root_i] = root_j;
                volume[root_j] += volume[root_i]; 
            }
        }
    };

    UnionFind uf(n);
    float max_radius = 0.0f;
    double sum_radius = 0.0;

    // Parallelize Stats Collection
    #pragma omp parallel for reduction(max:max_radius) reduction(+:sum_radius)
    for (int i = 0; i < n; ++i) {
        float r = sphere_table(i, 3);
        uf.volume[i] = r * r * r; 
        sum_radius += r;
        if (r > max_radius) max_radius = r;
    }

    const float voxel_diag = 1.733f; 
    float avg_radius = static_cast<float>(sum_radius / n);
    float cell_size = std::max(1.0f, (2.0f * avg_radius) + voxel_diag);

    Eigen::Vector3f min_b = sphere_table.leftCols<3>().colwise().minCoeff();
    Eigen::Vector3f max_b = sphere_table.leftCols<3>().colwise().maxCoeff();
    
    int grid_nx = std::max(1, static_cast<int>(std::ceil((max_b.x() - min_b.x()) / cell_size)) + 1);
    int grid_ny = std::max(1, static_cast<int>(std::ceil((max_b.y() - min_b.y()) / cell_size)) + 1);
    int grid_nz = std::max(1, static_cast<int>(std::ceil((max_b.z() - min_b.z()) / cell_size)) + 1);
    
    long long total_cells_ll = static_cast<long long>(grid_nx) * grid_ny * grid_nz;
    const long long MAX_SAFE_CELLS = 16777216; 
    
    if (total_cells_ll > MAX_SAFE_CELLS) {
        throw std::runtime_error("[Dev-Multisphere] Spatial hash grid allocation exceeded safety threshold.");
    }

    int total_cells = static_cast<int>(total_cells_ll);
    std::vector<int> head(total_cells, -1);
    std::vector<int> next(n, -1);

    auto get_cell_coords = [&](float x, float y, float z, int& cx, int& cy, int& cz) {
        cx = std::max(0, std::min(grid_nx - 1, static_cast<int>((x - min_b.x()) / cell_size)));
        cy = std::max(0, std::min(grid_ny - 1, static_cast<int>((y - min_b.y()) / cell_size)));
        cz = std::max(0, std::min(grid_nz - 1, static_cast<int>((z - min_b.z()) / cell_size)));
    };

    auto get_1d_idx = [&](int cx, int cy, int cz) -> int {
        return cx * (grid_ny * grid_nz) + cy * grid_nz + cz;
    };

    // Serial grid population (ultra-fast, lock-free continuous memory)
    for (int i = 0; i < n; ++i) {
        int cx, cy, cz;
        get_cell_coords(sphere_table(i, 0), sphere_table(i, 1), sphere_table(i, 2), cx, cy, cz);
        int cell_idx = get_1d_idx(cx, cy, cz);
        next[i] = head[cell_idx];
        head[cell_idx] = i;
    }

    // --- PHASE 1: PARALLEL SEARCH ---
    std::vector<std::pair<int, int>> global_edges;

    #pragma omp parallel
    {
        std::vector<std::pair<int, int>> local_edges;
        
        #pragma omp for schedule(dynamic, 64)
        for (int i = 0; i < n; ++i) {
            float xi = sphere_table(i, 0), yi = sphere_table(i, 1), zi = sphere_table(i, 2), ri = sphere_table(i, 3);
            int cx, cy, cz;
            get_cell_coords(xi, yi, zi, cx, cy, cz);
            int search_extent = static_cast<int>(std::ceil((ri + max_radius + voxel_diag) / cell_size));

            for (int dx = -search_extent; dx <= search_extent; ++dx) {
                int nx = cx + dx;
                if (nx < 0 || nx >= grid_nx) continue;
                for (int dy = -search_extent; dy <= search_extent; ++dy) {
                    int ny = cy + dy;
                    if (ny < 0 || ny >= grid_ny) continue;
                    for (int dz = -search_extent; dz <= search_extent; ++dz) {
                        int nz = cz + dz;
                        if (nz < 0 || nz >= grid_nz) continue;

                        int neighbor_cell = get_1d_idx(nx, ny, nz);
                        int j = head[neighbor_cell];
                        
                        while (j != -1) {
                            if (i < j) { 
                                float rj = sphere_table(j, 3);
                                float dist_sq = (xi - sphere_table(j, 0)) * (xi - sphere_table(j, 0)) + 
                                                (yi - sphere_table(j, 1)) * (yi - sphere_table(j, 1)) + 
                                                (zi - sphere_table(j, 2)) * (zi - sphere_table(j, 2));
                                float rad_sum = ri + rj + voxel_diag;
                                
                                if (dist_sq <= rad_sum * rad_sum) {
                                    local_edges.emplace_back(i, j);
                                }
                            }
                            j = next[j];
                        }
                    }
                }
            }
        }
        
        #pragma omp critical
        {
            global_edges.insert(global_edges.end(), local_edges.begin(), local_edges.end());
        }
    }

    // --- PHASE 2: SERIAL UNION ---
    for (const auto& edge : global_edges) {
        uf.unite(edge.first, edge.second);
    }

    // Find root with maximum accumulated volume proxy
    int max_root = -1;
    float max_vol = -1.0f;
    for (int i = 0; i < n; ++i) {
        int root = uf.find(i);
        if (uf.volume[root] > max_vol) {
            max_vol = uf.volume[root];
            max_root = root;
        }
    }

    // Extract largest network
    std::vector<Eigen::Vector4f> kept_spheres;
    kept_spheres.reserve(n);
    for (int i = 0; i < n; ++i) {
        if (uf.find(i) == max_root) {
            kept_spheres.push_back(sphere_table.row(i));
        }
    }

    Eigen::MatrixX4f result(kept_spheres.size(), 4);
    
    #pragma omp parallel for if(kept_spheres.size() > 1000)
    for (int i = 0; i < static_cast<int>(kept_spheres.size()); ++i) {
        result.row(i) = kept_spheres[i];
    }

    return result;
}


/**
 * @brief Appends peaks to the sphere table, respecting max_spheres limit.
 * @param table Sphere table (modified in place).
 * @param peaks Peaks to append.
 * @param max_spheres Optional maximum number of spheres.
 */
inline void append_sphere_table(
    Eigen::MatrixX4f& table,
    const Eigen::MatrixX4f& peaks,
    int max_spheres = 0)
{
    if (peaks.rows() == 0) return;
    if (max_spheres == 0 ||
        (table.rows() + peaks.rows() < static_cast<size_t>(max_spheres))) {
        table.conservativeResize(table.rows() + peaks.rows(), 4);
        for (int i = 0; i < peaks.rows(); ++i) {
            table.row(table.rows() - peaks.rows() + i) = peaks.row(i);
        }
        return;
    } else {
        int space_left = max_spheres - table.rows();
        if (space_left <= 0) return;
        table.conservativeResize(table.rows() + space_left, 4);
        for (int i = 0; i < space_left; ++i) {
            table.row(table.rows() - space_left + i) = peaks.row(i);
        }
    }
}


/**
 * @brief Summed field F = D + weight * max(0, D - R), evaluated where it is read instead of stored.
 *
 * D is the distance field of the target, R that of the reconstruction. With weight 0, F is exactly D.
 */
struct SummedField {
    const float* original;        ///< D, distance field of the target
    const float* reconstruction;  ///< R, distance field of the reconstruction
    float weight;                 ///< weight of the residual max(0, D - R)

    /// Field value at linear index i.
    float operator[](long i) const {
        float residual = original[i] - reconstruction[i];
        residual = (residual > 0.0f) ? residual : 0.0f;
        return original[i] + residual * weight;
    }
};

/**
 * @brief Iterative multisphere solver.
 *
 * Each iteration places the peaks of the summed field as spheres, until the precision target is reached or no
 * peaks are found `persistence` times in a row (each time with a larger residual weight). The solver owns its
 * grids and passes change information between the stages instead of rescanning the whole grid:
 *   - rasterizing new spheres records every voxel that changes from 0 to 1: it marks the voxel's line for the
 *     incremental distance transform and updates the Dice counters;
 *   - the distance field of the reconstruction is updated incrementally (IncrementalEDT), which reports the
 *     blocks in which it changed;
 *   - the summed field is computed where it is read (SummedField);
 *   - while the weight is unchanged, the peak search re-tests only the blocks around changed distances.
 * The result is bit-identical to recomputing every stage on the full grid.
 */
class MultisphereSolver {
public:
    /**
     * @brief Prepares a solver.
     * @param original_distance Distance field of the target.
     * @param target Target voxel grid (> 0 is inside).
     * @param recon_mask Reconstruction mask; receives the reconstruction (usually empty at the start).
     * @param config Solver settings (search_window, min_radius_vox, min_center_distance_rel, radius_offset_vox,
     *        precision_target, persistence, max_spheres, initial_sphere_table, show_progress).
     * @throws std::invalid_argument if the grid shapes differ.
     */
    MultisphereSolver(const VoxelGrid<float>& original_distance,
                      const VoxelGrid<uint8_t>& target,
                      VoxelGrid<uint8_t>& recon_mask,
                      const MultisphereConfig& config)
        : original_distance_(original_distance),
          target_(target),
          recon_mask_(recon_mask),
          config_(config),
          overlap_mask_(target.nx(), target.ny(), target.nz(), target.voxel_size, target.origin),
          recon_distance_(target.nx(), target.ny(), target.nz(), target.voxel_size, target.origin),
          recon_edt_(static_cast<int>(target.nx()), static_cast<int>(target.ny()), static_cast<int>(target.nz()),
                     std::max(8, config.search_window))  // blocks at least as wide as the search window
    {
        if (recon_mask.shape != target.shape || original_distance.shape != target.shape) {
            throw std::invalid_argument("MultisphereSolver: grid shapes must match.");
        }
        size_t target_volume = 0, recon_volume = 0, intersection = 0;
        #pragma omp parallel for reduction(+:target_volume, recon_volume, intersection)
        for (size_t i = 0; i < target.data.size(); ++i) {
            const bool t = target.data[i] > 0;
            const bool r = recon_mask.data[i] > 0;
            target_volume += t;
            recon_volume += r;
            intersection += (t && r);
        }
        target_volume_ = target_volume;
        recon_volume_ = recon_volume;
        intersection_ = intersection;
        test_blocks_.assign(recon_edt_.changed_blocks().size(), 0);
    }

    /**
     * @brief Runs the solver.
     * @return Sphere table (Nx4: x, y, z, radius in voxel units), starting with config.initial_sphere_table.
     *         recon_mask holds the union of the spheres.
     */
    Eigen::MatrixX4f run() {
        Eigen::MatrixX4f sphere_table = config_.initial_sphere_table;
        int iter = 1;                // attempts since the last iteration that found peaks
        int prev_count = 0;          // number of spheres before the last accepted batch
        float weight_factor = 1.0f;  // residual weight, increased after each attempt without peaks
        bool peaks_found = false;    // the last iteration accepted spheres

        if (sphere_table.rows() > 0) {
            rasterize_reconstruction(sphere_table);
            prev_count = sphere_table.rows();
            rasterize_overlap(sphere_table);
            peaks_found = true;
        }

        while (iter <= config_.persistence) {
            if (config_.max_spheres > 0 && sphere_table.rows() >= config_.max_spheres) {
                if (config_.show_progress) std::cout << "Reached maximum number of spheres. Terminating." << std::endl;
                break;
            }

            // A. Check the precision and update the fields after a new batch of spheres
            float field_weight = 0.0f;  // no spheres yet: the summed field is the original distance field
            bool distance_updated = false;
            if (sphere_table.rows() > 0) {
                if (peaks_found) {
                    iter = 1;
                    const float precision = static_cast<float>(dice());
                    if (config_.show_progress) {
                        std::cout << " -Weight = " << weight_factor << " -Total spheres = " << sphere_table.rows()
                                  << " -Precision = " << precision << std::endl;
                    }
                    if (config_.precision_target < precision) {
                        if (config_.show_progress) std::cout << "Reached target precision. Terminating." << std::endl;
                        break;
                    }
                    rasterize_overlap(sphere_table.block(prev_count, 0, sphere_table.rows() - prev_count, 4));
                    recon_edt_.update(recon_mask_, recon_distance_);
                    distance_updated = true;
                }
                field_weight = weight_factor;
            }
            const SummedField field{original_distance_.data.data(), recon_distance_.data.data(), field_weight};

            // B. Find and filter the peaks
            Eigen::MatrixX4f peaks = find_peaks(field, distance_updated);
            peaks = filter_and_shift_peaks(peaks, sphere_table, field, static_cast<int>(target_.nx()),
                                           static_cast<int>(target_.ny()), static_cast<int>(target_.nz()), config_);

            // C. Accept the peaks as spheres, or retry with a larger weight
            if (peaks.rows() == 0) {
                peaks_found = false;
                ++iter;
                ++weight_factor;
                continue;
            }
            peaks_found = true;
            prev_count = sphere_table.rows();
            append_sphere_table(sphere_table, peaks, config_.max_spheres);
            rasterize_reconstruction(sphere_table.block(prev_count, 0, sphere_table.rows() - prev_count, 4));
        }
        return sphere_table;
    }

private:
    const VoxelGrid<float>& original_distance_;  ///< distance field of the target
    const VoxelGrid<uint8_t>& target_;           ///< target voxel grid
    VoxelGrid<uint8_t>& recon_mask_;             ///< union of the accepted spheres
    MultisphereConfig config_;

    VoxelGrid<uint8_t> overlap_mask_;  ///< voxels near accepted spheres, excluded from the peak search
    VoxelGrid<float> recon_distance_;  ///< distance field of the reconstruction, updated by recon_edt_
    IncrementalEDT recon_edt_;

    size_t target_volume_ = 0;  ///< target voxels
    size_t recon_volume_ = 0;   ///< reconstruction voxels
    size_t intersection_ = 0;   ///< voxels in both

    std::vector<Eigen::Vector4f> peak_cache_;  ///< peaks of the last search
    bool peak_cache_valid_ = false;            ///< a search has run
    float peak_cache_weight_ = 0.0f;           ///< summed-field weight of the last search
    std::vector<uint8_t> test_blocks_;         ///< blocks to re-test in an incremental search

    /**
     * @brief Peaks of the summed field.
     *
     * While the weight is unchanged, the field changes only where the reconstruction distance changed, so only
     * those blocks and their neighbours are re-tested (peak_local_max_3d_incremental). A full search runs on the
     * first call, after a weight change, after the first distance transform, or when most blocks changed.
     *
     * @param field Summed field of this iteration.
     * @param distance_updated The reconstruction distance was updated since the last search.
     * @return Peaks (Nx4: x, y, z, radius), as from peak_local_max_3d().
     */
    Eigen::MatrixX4f find_peaks(const SummedField& field, bool distance_updated) {
        bool full = !peak_cache_valid_ || field.weight != peak_cache_weight_;
        if (!full) {
            if (distance_updated) {
                full = recon_edt_.all_changed() || !mark_test_blocks();
            } else {
                std::fill(test_blocks_.begin(), test_blocks_.end(), 0);  // field unchanged: only the mask grew
            }
        }
        peak_cache_valid_ = true;
        peak_cache_weight_ = field.weight;

        if (full) {
            Eigen::MatrixX4f peaks = peak_local_max_3d(field, original_distance_, overlap_mask_, config_);
            peak_cache_.resize(peaks.rows());
            for (int i = 0; i < peaks.rows(); ++i) peak_cache_[i] = peaks.row(i).transpose();
            return peaks;
        }
        return peak_local_max_3d_incremental(field, original_distance_, overlap_mask_, config_,
                                             test_blocks_, recon_edt_.block_size(), peak_cache_);
    }

    /**
     * @brief Marks the blocks to re-test: every block with a changed distance and its 26 neighbours. Blocks are
     *        at least as wide as the search window, so this covers every window that contains a change.
     * @return False if more than half of the blocks are marked (a full search is then cheaper).
     */
    bool mark_test_blocks() {
        const int nbx = recon_edt_.block_grid()[0], nby = recon_edt_.block_grid()[1], nbz = recon_edt_.block_grid()[2];
        const std::vector<uint8_t>& changed = recon_edt_.changed_blocks();
        std::fill(test_blocks_.begin(), test_blocks_.end(), 0);
        for (int bx = 0; bx < nbx; ++bx) {
            for (int by = 0; by < nby; ++by) {
                for (int bz = 0; bz < nbz; ++bz) {
                    if (!changed[(static_cast<long>(bx) * nby + by) * nbz + bz]) continue;
                    for (int x = std::max(0, bx - 1); x <= std::min(nbx - 1, bx + 1); ++x)
                        for (int y = std::max(0, by - 1); y <= std::min(nby - 1, by + 1); ++y)
                            for (int z = std::max(0, bz - 1); z <= std::min(nbz - 1, bz + 1); ++z)
                                test_blocks_[(static_cast<long>(x) * nby + y) * nbz + z] = 1;
                }
            }
        }
        size_t marked = 0;
        for (uint8_t flag : test_blocks_) marked += flag;
        return marked * 2 <= test_blocks_.size();
    }

    /**
     * @brief Renders spheres into the reconstruction mask (the same voxels as spheres_to_grid()).
     *
     * Every voxel that changes from 0 to 1 marks its x-line for the incremental distance transform and updates
     * the Dice counters. Spheres overlap and are rendered in parallel, so each voxel is claimed atomically and
     * counted once.
     *
     * @param spheres Spheres (Nx4: x, y, z, radius in voxel units).
     */
    void rasterize_reconstruction(const Eigen::MatrixX4f& spheres) {
        if (spheres.rows() == 0) return;
        const int n_x = static_cast<int>(recon_mask_.nx());
        const int n_y = static_cast<int>(recon_mask_.ny());
        const int n_z = static_cast<int>(recon_mask_.nz());
        const long stride_x = static_cast<long>(n_y) * n_z;
        uint8_t* mask = recon_mask_.data.data();
        const uint8_t* target = target_.data.data();

        size_t added = 0, added_inside = 0;
        #pragma omp parallel for reduction(+:added, added_inside)
        for (int i = 0; i < spheres.rows(); ++i) {
            const float cx = spheres(i, 0), cy = spheres(i, 1), cz = spheres(i, 2), radius = spheres(i, 3);
            if (radius <= 0) continue;
            const float r_sq = radius * radius;
            const int min_x = std::max(0, static_cast<int>(std::floor(cx - radius)));
            const int max_x = std::min(n_x - 1, static_cast<int>(std::ceil(cx + radius)));
            const int min_y = std::max(0, static_cast<int>(std::floor(cy - radius)));
            const int max_y = std::min(n_y - 1, static_cast<int>(std::ceil(cy + radius)));
            const int min_z = std::max(0, static_cast<int>(std::floor(cz - radius)));
            const int max_z = std::min(n_z - 1, static_cast<int>(std::ceil(cz + radius)));
            for (int x = min_x; x <= max_x; ++x) {
                const float dx = x - cx;
                const float dx2 = dx * dx;
                if (dx2 > r_sq) continue;
                for (int y = min_y; y <= max_y; ++y) {
                    const float dy = y - cy;
                    const float dxy2 = dx2 + dy * dy;
                    if (dxy2 > r_sq) continue;
                    for (int z = min_z; z <= max_z; ++z) {
                        const float dz = z - cz;
                        if (dxy2 + dz * dz > r_sq) continue;
                        const long idx = x * stride_x + static_cast<long>(y) * n_z + z;
                        if (__atomic_load_n(&mask[idx], __ATOMIC_RELAXED) != 0) continue;
                        if (__atomic_exchange_n(&mask[idx], static_cast<uint8_t>(1), __ATOMIC_RELAXED) != 0) continue;
                        ++added;
                        added_inside += target[idx] > 0;
                        recon_edt_.mark_line(y, z);
                    }
                }
            }
        }
        recon_volume_ += added;
        intersection_ += added_inside;
    }

    /**
     * @brief Renders spheres into the overlap mask, with radius sqrt(r) * min_center_distance_rel.
     * @param spheres Spheres (Nx4: x, y, z, radius in voxel units).
     */
    void rasterize_overlap(const Eigen::MatrixX4f& spheres) {
        spheres_to_grid<uint8_t>(overlap_mask_, spheres, 1, config_.min_center_distance_rel);
    }

    /// Dice coefficient of reconstruction and target from the counters (same formula as compute_voxel_precision()).
    double dice() const {
        if (target_volume_ == 0 && recon_volume_ == 0) return 1.0;
        if (target_volume_ == 0 || recon_volume_ == 0) return 0.0;
        return (2.0 * static_cast<double>(intersection_)) / static_cast<double>(target_volume_ + recon_volume_);
    }
};

/**
 * @brief Iterative solver that places spheres on the peaks of the distance and residual fields.
 * @param original_distance Distance field of the target.
 * @param voxel_grid Target voxel grid (> 0 is inside).
 * @param recon_mask Reconstruction mask; receives the union of the spheres.
 * @param config Solver settings (see MultisphereSolver).
 * @return Sphere table (Nx4: x, y, z, radius in voxel units).
 */
inline Eigen::MatrixX4f compute_sphere_table(
    const VoxelGrid<float>& original_distance,
    const VoxelGrid<uint8_t>& voxel_grid,
    VoxelGrid<uint8_t>& recon_mask,
    MultisphereConfig config)
{
    return MultisphereSolver(original_distance, voxel_grid, recon_mask, config).run();
}

} // namespace GEMSS

#endif // GEMSS_RECONSTRUCTION_HELPERS_HPP
