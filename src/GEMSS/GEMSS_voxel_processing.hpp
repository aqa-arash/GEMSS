#ifndef GEMSS_VOXEL_PROCESSING_HPP
#define GEMSS_VOXEL_PROCESSING_HPP


/**
 * @file multisphere_voxel_processing.hpp
 * @brief Voxel grid processing utilities for GEMSS.
 *
 * Provides kernel application, peak shifting, precision computation, sphere rendering, and mesh conversion.
 *
 * @author Arash Moradian
 * @date 2026-03-09
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <Eigen/Core>
#ifdef HAVE_OPENMP
    #include <omp.h>
#endif

#include "GEMSS_datatypes.hpp"

namespace GEMSS {

/**
 * @brief Kernel operation modes for voxel grid processing.
 */
enum class KernelMode { ADD, SUBTRACT, ZERO };

/**
 * @brief Apply a kernel to a voxel grid at a specified center.
 * @tparam T Voxel grid data type.
 * @tparam K Kernel data type.
 * @param grid Target voxel grid.
 * @param center Center coordinate.
 * @param kernel Kernel grid.
 * @param mode Operation mode (add, subtract, zero).
 * @param scale Scaling factor for kernel values.
 */
template<typename T, typename K>
void apply_kernel_to_grid(
    VoxelGrid<T>& grid,
    const Eigen::Vector3f& center,
    const VoxelGrid<K>& kernel,
    KernelMode mode = KernelMode::ADD,
    float scale = 1.0f)
{
    const int dim = 3;
    Eigen::Vector3i k_shape(kernel.nx(), kernel.ny(), kernel.nz());
    Eigen::Vector3i v_shape(grid.nx(), grid.ny(), grid.nz());

    for (int i = 0; i < dim; ++i) {
        double floor_val;
        double frac = std::modf(center[i], &floor_val);
        if (k_shape[i] % 2 == 0) {
            if (std::abs(frac - 0.5) > 1e-7) throw std::invalid_argument("Even kernel requires .5 center.");
        } else {
            if (std::abs(frac) > 1e-7) throw std::invalid_argument("Odd kernel requires integer center.");
        }
    }

    Eigen::Vector3f half = (k_shape.cast<float>().array() - 1.0f) / 2.0f;
    Eigen::Vector3i llf_k = (center - half).array().round().cast<int>();
    Eigen::Vector3i trr_k = (center + half).array().round().cast<int>();

    Eigen::Vector3i start_v = llf_k.cwiseMax(0);
    Eigen::Vector3i end_v = trr_k.cwiseMin((v_shape.array() - 1).matrix());
    Eigen::Vector3i k_offset = start_v - llf_k;

    #pragma omp parallel for collapse(2) if((end_v[0] - start_v[0]) > 32)
    for (int x = start_v[0]; x <= end_v[0]; ++x) {
        for (int y = start_v[1]; y <= end_v[1]; ++y) {
            int vk_x = x - start_v[0] + k_offset[0];
            int vk_y = y - start_v[1] + k_offset[1];
            for (int z = start_v[2]; z <= end_v[2]; ++z) {
                int vk_z = z - start_v[2] + k_offset[2];
                T& vol_val = grid(x, y, z);
                const T& ker_val = kernel.data[vk_x * (kernel.ny() * kernel.nz()) + vk_y * kernel.nz() + vk_z];
                switch (mode) {
                    case KernelMode::ADD:
                        vol_val += static_cast<T>(ker_val * scale);
                        break;
                    case KernelMode::SUBTRACT:
                        vol_val -= static_cast<T>(ker_val * scale);
                        break;
                    case KernelMode::ZERO:
                        if (ker_val != 0) vol_val = 0;
                        break;
                }
            }
        }
    }
}

/**
 * @brief Helper function to check axis direction for peak shifting.
 * @param field Voxel grid field.
 * @param cx Center x.
 * @param cy Center y.
 * @param cz Center z.
 * @param dx Direction x.
 * @param dy Direction y.
 * @param dz Direction z.
 * @param current_val Current value.
 * @param epsilon Tolerance.
 * @return Shift value (1.0, -1.0, or 0.0).
 */
inline float check_axis_shift(
    const VoxelGrid<float>& field,
    int cx, int cy, int cz,
    int dx, int dy, int dz,
    float current_val,
    float epsilon)
{
    int nx = field.nx();
    int ny = field.ny();
    int nz = field.nz();

    int n1_x = cx + dx;
    int n1_y = cy + dy;
    int n1_z = cz + dz;
    if (n1_x < 0 || n1_x >= nx || n1_y < 0 || n1_y >= ny || n1_z < 0 || n1_z >= nz) return 0.0f;
    float val_n1 = field(n1_x, n1_y, n1_z);
    if (std::abs(val_n1 - current_val) > epsilon) return 0.0f;

    int n2_x = cx + (dx * 2);
    int n2_y = cy + (dy * 2);
    int n2_z = cz + (dz * 2);
    if (n2_x < 0 || n2_x >= nx || n2_y < 0 || n2_y >= ny || n2_z < 0 || n2_z >= nz) return (float)dx;
    float val_n2 = field(n2_x, n2_y, n2_z);
    if (val_n2 < current_val - epsilon) {
        return (float)dx;
    }
    return 0.0f;
}

/**
 * @brief Determines sub-voxel shift for local maxima.
 * @param field Voxel grid field.
 * @param idx Voxel index.
 * @param epsilon Tolerance.
 * @return Sub-voxel shift vector.
 */
/**
 * @brief Determines sub-voxel shift for local maxima utilizing a zero-bounds-check fast path.
 */
inline Eigen::Vector3f shift_voxel_center(
    const VoxelGrid<float>& field,
    const Eigen::Vector3i& idx,
    float epsilon = 1e-5f)
{
    int x = idx.x(), y = idx.y(), z = idx.z();
    int nx = static_cast<int>(field.nx());
    int ny = static_cast<int>(field.ny());
    int nz = static_cast<int>(field.nz());

    // FAST-PATH: Strictly interior voxels (avoids 36 conditional bounds checks)
    if (x >= 2 && x < nx - 2 && y >= 2 && y < ny - 2 && z >= 2 && z < nz - 2) {
        const float* data_ptr = field.data.data();
        size_t stride_y = nz;
        size_t stride_x = ny * nz;
        size_t base_idx = x * stride_x + y * stride_y + z;

        float val = data_ptr[base_idx];
        float sx = 0.0f, sy = 0.0f, sz = 0.0f;

        // X-axis check (Prefetched contiguous strides)
        if (std::abs(data_ptr[base_idx + stride_x] - val) <= epsilon && 
            data_ptr[base_idx + 2 * stride_x] < val - epsilon) sx = 1.0f;
        else if (std::abs(data_ptr[base_idx - stride_x] - val) <= epsilon && 
                 data_ptr[base_idx - 2 * stride_x] < val - epsilon) sx = -1.0f;

        // Y-axis check
        if (std::abs(data_ptr[base_idx + stride_y] - val) <= epsilon && 
            data_ptr[base_idx + 2 * stride_y] < val - epsilon) sy = 1.0f;
        else if (std::abs(data_ptr[base_idx - stride_y] - val) <= epsilon && 
                 data_ptr[base_idx - 2 * stride_y] < val - epsilon) sy = -1.0f;

        // Z-axis check (Perfect cache locality)
        if (std::abs(data_ptr[base_idx + 1] - val) <= epsilon && 
            data_ptr[base_idx + 2] < val - epsilon) sz = 1.0f;
        else if (std::abs(data_ptr[base_idx - 1] - val) <= epsilon && 
                 data_ptr[base_idx - 2] < val - epsilon) sz = -1.0f;

        return Eigen::Vector3f(sx, sy, sz);
    }

    // SLOW-PATH: Boundary voxels fallback
    float val = field(x, y, z);
    float sx = check_axis_shift(field, x, y, z, 1, 0, 0, val, epsilon);
    if (sx == 0.0f) sx = check_axis_shift(field, x, y, z, -1, 0, 0, val, epsilon);

    float sy = check_axis_shift(field, x, y, z, 0, 1, 0, val, epsilon);
    if (sy == 0.0f) sy = check_axis_shift(field, x, y, z, 0, -1, 0, val, epsilon);

    float sz = check_axis_shift(field, x, y, z, 0, 0, 1, val, epsilon);
    if (sz == 0.0f) sz = check_axis_shift(field, x, y, z, 0, 0, -1, val, epsilon);

    return Eigen::Vector3f(sx, sy, sz);
}

/**
 * @brief Computes reconstruction precision based on mismatch.
 * @tparam T Target voxel grid data type.
 * @tparam U Reconstruction voxel grid data type.
 * @param target Target voxel grid.
 * @param reconstruction Reconstruction voxel grid.
 * @return Precision value [0, 1].
 */
template <typename T, typename U>
double compute_voxel_precision(const VoxelGrid<T>& target,
                               const VoxelGrid<U>& reconstruction)
{
    if (target.nx() != reconstruction.nx() || target.ny() != reconstruction.ny() || target.nz() != reconstruction.nz()) {
        throw std::invalid_argument("Shapes must match.");
    }

    size_t total_target = 0;
    size_t mismatches = 0;
    const size_t n = target.data.size();

    #pragma omp parallel for reduction(+:total_target, mismatches)
    for (size_t i = 0; i < n; ++i) {
        bool t = (target.data[i] > static_cast<T>(0));
        bool r = (reconstruction.data[i] > static_cast<U>(0));
        if (t) total_target++;
        if (t != r) mismatches++;
    }

    if (total_target == 0) {
        size_t total_recon = 0;
        #pragma omp parallel for reduction(+:total_recon)
        for (size_t i = 0; i < n; ++i) if (reconstruction.data[i]) total_recon++;
        return (total_recon == 0) ? 1.0 : 0.0;
    }

    float mismatch_fraction = static_cast<float>(mismatches) / static_cast<float>(total_target);
    return std::clamp(1.0f - mismatch_fraction, 0.0f, 1.0f);
}

/**
 * @brief Renders spheres into a voxel grid.
 * @tparam T Voxel grid data type.
 * @param grid Target voxel grid.
 * @param sphere_table Sphere table (Nx4 matrix: x, y, z, diameter_vox).
 * @param fill_value Value to fill inside spheres.
 * @param coefficient Coefficient for adjusting sphere sizes.
 */
template <typename T>
void spheres_to_grid(VoxelGrid<T>& grid,
    const Eigen::MatrixX4f& sphere_table,
    T fill_value = static_cast<T>(1), float coefficient = 1.0f)
{
    if (sphere_table.rows() == 0) return;
    #pragma omp parallel for
    for (int i = 0; i < sphere_table.rows(); ++i) {
        float cx = sphere_table(i, 0);
        float cy = sphere_table(i, 1);
        float cz = sphere_table(i, 2);
        float radius_vox = sphere_table(i, 3) * coefficient;
        if (radius_vox <= 0) continue;
        grid.sphere_kernel(cx, cy, cz, radius_vox, fill_value);
    }
}

/**
 * @brief Converts a voxel grid into a blocky mesh by generating faces for exposed voxels.
 * @tparam T Voxel grid data type.
 * @param grid Input voxel grid.
 * @param threshold Threshold for voxel occupancy.
 * @return FastMesh mesh structure.
 */
template <typename T>
FastMesh grid_to_mesh(const VoxelGrid<T>& grid, T threshold = static_cast<T>(0)) {
    std::vector<Eigen::Vector3f> out_verts;
    std::vector<Eigen::Vector3i> out_tris;

    out_verts.reserve(grid.data.size() / 5);
    out_tris.reserve(grid.data.size() / 5);

    int nx = grid.nx();
    int ny = grid.ny();
    int nz = grid.nz();
    float vs = grid.voxel_size;

    const int dx[6] = {1, -1, 0, 0, 0, 0};
    const int dy[6] = {0, 0, 1, -1, 0, 0};
    const int dz[6] = {0, 0, 0, 0, 1, -1};

    const float face_verts[6][4][3] = {
        {{1,0,0}, {1,1,0}, {1,1,1}, {1,0,1}}, // +x
        {{0,0,0}, {0,0,1}, {0,1,1}, {0,1,0}}, // -x
        {{0,1,0}, {0,1,1}, {1,1,1}, {1,1,0}}, // +y
        {{0,0,0}, {1,0,0}, {1,0,1}, {0,0,1}}, // -y
        {{0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}}, // +z
        {{0,0,0}, {0,1,0}, {1,1,0}, {1,0,0}}  // -z
    };

    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int z = 0; z < nz; ++z) {
                if (grid(x, y, z) <= threshold) continue;
                for (int f = 0; f < 6; ++f) {
                    int nx_idx = x + dx[f];
                    int ny_idx = y + dy[f];
                    int nz_idx = z + dz[f];
                    bool draw_face = false;
                    if (nx_idx < 0 || nx_idx >= nx ||
                        ny_idx < 0 || ny_idx >= ny ||
                        nz_idx < 0 || nz_idx >= nz) {
                        draw_face = true;
                    }
                    else if (grid(nx_idx, ny_idx, nz_idx) <= threshold) {
                        draw_face = true;
                    }
                    if (draw_face) {
                        int base_idx = static_cast<int>(out_verts.size());
                        for (int v = 0; v < 4; ++v) {
                            Eigen::Vector3f pos;
                            pos.x() = grid.origin.x() + (x + face_verts[f][v][0]) * vs;
                            pos.y() = grid.origin.y() + (y + face_verts[f][v][1]) * vs;
                            pos.z() = grid.origin.z() + (z + face_verts[f][v][2]) * vs;
                            out_verts.push_back(pos);
                        }
                        out_tris.push_back({base_idx, base_idx + 1, base_idx + 2});
                        out_tris.push_back({base_idx, base_idx + 2, base_idx + 3});
                    }
                }
            }
        }
    }

    FastMesh mesh;
    mesh.vertices.resize(out_verts.size(), 3);
    mesh.triangles.resize(out_tris.size(), 3);

    #pragma omp parallel for
    for (size_t i = 0; i < out_verts.size(); ++i) mesh.vertices.row(i) = out_verts[i];
    #pragma omp parallel for
    for (size_t i = 0; i < out_tris.size(); ++i) mesh.triangles.row(i) = out_tris[i];

    return mesh;
}

} // namespace MSS

#endif // GEMSS_VOXEL_PROCESSING_HPP
