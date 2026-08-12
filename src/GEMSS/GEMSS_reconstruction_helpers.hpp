
/**
 * @file multisphere_reconstruction_helpers.hpp
 * @brief Helper functions for GEMSS reconstruction algorithms.
 *
 * Provides peak detection, filtering, sphere table management, and residual field computation.
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

namespace GEMSS {

/**
 * @brief Detects local maxima in a 3D volume.
 * @param field Distance field to search for peaks.
 * @param original_distance Original distance field for radius values.
 * @param min_distance Minimum distance between peaks.
 * @param min_radius_vox Minimum radius threshold.
 * @return Nx4 matrix (x, y, z, radius) of detected peaks.
 */
inline Eigen::MatrixX4f peak_local_max_3d(
    const VoxelGrid<float>& field,
    const VoxelGrid<float>& original_distance,
    const VoxelGrid<uint8_t>& overlap_mask,
    MultisphereConfig config)
{
    std::vector<Eigen::Vector4f> peaks;
    int nx = (int)field.nx();
    int ny = (int)field.ny();
    int nz = (int)field.nz();

    #pragma omp parallel
    {
        std::vector<Eigen::Vector4f> local_peaks;
        #pragma omp for collapse(2)
        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < ny; ++y) {
                for (int z = 0; z < nz; ++z) {
                    float val = field(x, y, z);
                    float radius = original_distance(x, y, z) ;

                    if (radius < config.min_radius_vox) continue;
                    if (val <= 1.0f) continue;
                    if (overlap_mask(x, y, z) > 0) continue; // Skip if within overlap region

                    bool is_max = true;
                    // Check local neighborhood
                    for (int dx = -config.search_window; is_max && dx <= config.search_window; ++dx) {
                        for (int dy = -config.search_window; is_max && dy <= config.search_window; ++dy) {
                            for (int dz = -config.search_window; is_max && dz <= config.search_window; ++dz) {
                                int nx_idx = x + dx, ny_idx = y + dy, nz_idx = z + dz;
                                if (nx_idx >= 0 && nx_idx < nx && ny_idx >= 0 && ny_idx < ny && nz_idx >= 0 && nz_idx < nz) {
                                    if (field(nx_idx, ny_idx, nz_idx) > val) is_max = false;
                                }
                            }
                        }
                    }
                    if (is_max) {
                        local_peaks.push_back(Eigen::Vector4f((float)x, (float)y, (float)z, radius + config.radius_offset_vox));
                    }
                }
            }
        }
        #pragma omp critical
        peaks.insert(peaks.end(), local_peaks.begin(), local_peaks.end());
    }

    // Sort based on radius (descending), then x,y,z (ascending)
    std::sort(peaks.begin(), peaks.end(), [](const Eigen::Vector4f& a, const Eigen::Vector4f& b) {
        float val_a = a[3];
        float val_b = b[3];
        const float epsilon = 1e-5f;
        if (std::abs(val_a - val_b) > epsilon) {
            return val_a > val_b;
        }
        if (a[0] != b[0]) return a[0] < b[0];
        if (a[1] != b[1]) return a[1] < b[1];
        return a[2] < b[2];
    });

    Eigen::MatrixX4f result(peaks.size(), 4);
    for (size_t i = 0; i < peaks.size(); ++i) {
        result.row(i) = peaks[i];
    }
    return result;
}

/**
 * @brief Filters peaks based on minimum center distance and applies sub-voxel shifts.
 * @param peaks Candidate peaks (Nx4 matrix).
 * @param sphere_table History table of accepted spheres.
 * @param distance_field Distance field for sub-voxel shifting.
 * @param config Configuration parameters.
 * @return Filtered peaks as Nx4 matrix (x, y, z, radius).
 */
inline Eigen::MatrixX4f filter_and_shift_peaks(
    const Eigen::MatrixX4f& peaks,
    const Eigen::MatrixX4f& sphere_table,
    VoxelGrid<float>& distance_field,
    MultisphereConfig config)
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
            Eigen::Vector3f direction = shift_voxel_center(distance_field, coord_int, 0.25f);
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
 * @brief Computes the residual distance field (original - spheres).
 * @param original_distance Original distance field.
 * @param spheres_distance Distance field from spheres.
 * @return Residual distance field.
 */
inline VoxelGrid<float> residual_distance_field(
    const VoxelGrid<float>& original_distance,
    const VoxelGrid<float>& spheres_distance)
{
    if (original_distance.nx() != spheres_distance.nx() ||
        original_distance.ny() != spheres_distance.ny() ||
        original_distance.nz() != spheres_distance.nz()) {
        throw std::invalid_argument("Shapes must match.");
    }
    VoxelGrid<float> residual(original_distance.nx(), original_distance.ny(), original_distance.nz(),
                              original_distance.voxel_size, original_distance.origin);
    #pragma omp parallel for
    for (size_t i = 0; i < residual.data.size(); ++i) {
        float val = original_distance.data[i] - spheres_distance.data[i];
        residual.data[i] = (val > 0.0f) ? val : 0.0f;
    }
    return residual;
}


/**
 * @brief Core iterative solver to place spheres based on distance fields and residuals.
 */
inline Eigen::MatrixX4f compute_sphere_table(
    const VoxelGrid<float>& original_distance,
    const VoxelGrid<uint8_t>& voxel_grid,
    VoxelGrid<uint8_t>& recon_mask,
    MultisphereConfig config
) {
    Eigen::MatrixX4f sphere_table = config.initial_sphere_table;

    float diff = static_cast<float>(config.search_window) - 1.0f;
    /* Based on off axis angle for 1 neighbor (stencil size 3) up to 12 is safe. */
    float stable_weight = 2.0f * diff * diff + 4.0f;

    // Solver State
    int max_iter = config.persistence, iter = 1, prev_count = 0;
    float weight_factor = 1.0f;
    bool peaks_found = false;

    VoxelGrid<float> running_field = original_distance;
    // =========================================================================
    // Pre-Loop Initialization
    // =========================================================================
    // If we start with an existing table, initialize the mask immediately.
    VoxelGrid<uint8_t> overlap_mask(voxel_grid.nx(), voxel_grid.ny(), voxel_grid.nz(), voxel_grid.voxel_size, voxel_grid.origin);
    if (sphere_table.rows() > 0) {
        spheres_to_grid<uint8_t>(recon_mask, sphere_table);
        prev_count = sphere_table.rows();
        spheres_to_grid<uint8_t>(overlap_mask, sphere_table, 1,  config.min_center_distance_rel);
        peaks_found = true;
    }

    if (peaks_found) {
        running_field.multiply(sphere_table.rows()* stable_weight +1);
        for (int i=0; i< sphere_table.rows(); ++i ){
            running_field.sphere_distance_kernel(sphere_table(i,0),sphere_table(i,1),sphere_table(i,2),sphere_table(i,3), stable_weight);
        }
    }
    
    // Temporary memory buffers scoped only to the solver
    VoxelGrid<float> residual(voxel_grid.nx(), voxel_grid.ny(), voxel_grid.nz(), voxel_grid.voxel_size, voxel_grid.origin);
    VoxelGrid<float> summed_field(voxel_grid.nx(), voxel_grid.ny(), voxel_grid.nz(), voxel_grid.voxel_size, voxel_grid.origin);
    

    while (iter <= max_iter) {
        if (config.max_spheres > 0  && sphere_table.rows() >= config.max_spheres) {
            if (config.show_progress) std::cout << "Reached maximum number of spheres. Terminating." << std::endl;
            break;
        }

        // A. Update Fields & Check Precision
        if (sphere_table.rows() > 0) {
            if (peaks_found) {
                iter = 1;
                float precision = compute_voxel_precision(voxel_grid, recon_mask);
                if (config.show_progress) {
                    std::cout << " -Weight = " << weight_factor << " -Total spheres = " << sphere_table.rows() << " -Precision = " << precision << std::endl;
                    }

                if (config.precision_target < precision) {
                    if (config.show_progress) std::cout << "Reached target precision. Terminating." << std::endl;
                    break;
                }
                spheres_to_grid<uint8_t>(overlap_mask, sphere_table.block(prev_count, 0, sphere_table.rows() - prev_count, 4), 1, config.min_center_distance_rel);
                //residual = residual_distance_field(original_distance, recon_mask.distance_transform());
                int added_spheres = sphere_table.rows() - prev_count;
                running_field.add(original_distance, added_spheres * stable_weight);
                for (int i = prev_count; i < sphere_table.rows();++i){
                    running_field.sphere_distance_kernel(sphere_table(i,0),sphere_table(i,1),sphere_table(i,2),sphere_table(i,3), stable_weight);
                }
                
            }
            #pragma omp parallel for
            for (size_t i = 0; i < summed_field.data.size(); ++i) {
                //summed_field.data[i] = original_distance.data[i] + residual.data[i] * weight_factor;

            // consider adding a running weight increase as per the old method 
            }
        } else {
            summed_field.data = original_distance.data;
        }

        // B. Peak Detection & Filtering
        //Eigen::MatrixX4f peaks = peak_local_max_3d(summed_field, original_distance, overlap_mask, config);
        Eigen::MatrixX4f peaks = peak_local_max_3d(running_field, original_distance, overlap_mask, config);
        peaks = filter_and_shift_peaks(peaks, sphere_table, summed_field,config);

        // C. Iteration Flow Control
        if (peaks.rows() == 0) {
            peaks_found = false;
            ++iter; 
            ++weight_factor;
            continue;
        }
        
        peaks_found = true;
        prev_count = sphere_table.rows();
        append_sphere_table(sphere_table, peaks, config.max_spheres);
        spheres_to_grid<uint8_t>(recon_mask, sphere_table.block(prev_count, 0, sphere_table.rows() - prev_count, 4));

    }

    return sphere_table;
}

} // namespace MSS

#endif // GEMSS_RECONSTRUCTION_HELPERS_HPP
