/**
 * @file multisphere_voxel_processing.hpp
 * @brief Voxel grid processing utilities for GEMSS.
 *
 * Provides kernel application, peak shifting, precision computation, sphere rendering, and mesh conversion.
 *
 * @author Arash Moradian
 * @date 2026-03-09
 */


#ifndef GEMSS_VOXEL_PROCESSING_HPP
#define GEMSS_VOXEL_PROCESSING_HPP


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
 * @brief Field accessor for a stored float grid, for the functions that accept any field (F[i], linear index,
 *        z fastest). See SummedField for a field that is computed on access.
 */
struct GridField {
    const float* data;  ///< grid values

    /// Value at linear index i.
    float operator[](long i) const { return data[i]; }
};

/**
 * @brief Sub-voxel shift of a local maximum.
 *
 * Along each axis, the shift is +1 (or -1) if the neighbour on that side has the same value (within epsilon)
 * and the next one is lower, i.e. the maximum is a plateau of two voxels. Interior voxels take a path without
 * bound checks; the boundary path follows check_axis_shift().
 *
 * @param F Field accessor, F[i] with i the linear index (z fastest).
 * @param nx, ny, nz Grid size.
 * @param idx Voxel of the maximum.
 * @param epsilon Tolerance for equal values.
 * @return Shift direction per axis (each -1, 0 or 1).
 */
template <typename Field>
inline Eigen::Vector3f shift_voxel_center(
    const Field& F,
    int nx, int ny, int nz,
    const Eigen::Vector3i& idx,
    float epsilon = 1e-5f)
{
    const int x = idx.x(), y = idx.y(), z = idx.z();
    const long stride_y = nz;
    const long stride_x = static_cast<long>(ny) * nz;
    const long base_idx = x * stride_x + y * stride_y + z;
    const float val = F[base_idx];

    // FAST-PATH: Strictly interior voxels (avoids 36 conditional bounds checks)
    if (x >= 2 && x < nx - 2 && y >= 2 && y < ny - 2 && z >= 2 && z < nz - 2) {
        float sx = 0.0f, sy = 0.0f, sz = 0.0f;

        // X-axis check
        if (std::abs(F[base_idx + stride_x] - val) <= epsilon &&
            F[base_idx + 2 * stride_x] < val - epsilon) sx = 1.0f;
        else if (std::abs(F[base_idx - stride_x] - val) <= epsilon &&
                 F[base_idx - 2 * stride_x] < val - epsilon) sx = -1.0f;

        // Y-axis check
        if (std::abs(F[base_idx + stride_y] - val) <= epsilon &&
            F[base_idx + 2 * stride_y] < val - epsilon) sy = 1.0f;
        else if (std::abs(F[base_idx - stride_y] - val) <= epsilon &&
                 F[base_idx - 2 * stride_y] < val - epsilon) sy = -1.0f;

        // Z-axis check
        if (std::abs(F[base_idx + 1] - val) <= epsilon &&
            F[base_idx + 2] < val - epsilon) sz = 1.0f;
        else if (std::abs(F[base_idx - 1] - val) <= epsilon &&
                 F[base_idx - 2] < val - epsilon) sz = -1.0f;

        return Eigen::Vector3f(sx, sy, sz);
    }

    // SLOW-PATH: Boundary voxels fallback (same rules as check_axis_shift)
    auto axis_shift = [&](int dx, int dy, int dz) -> float {
        const int n1_x = x + dx, n1_y = y + dy, n1_z = z + dz;
        if (n1_x < 0 || n1_x >= nx || n1_y < 0 || n1_y >= ny || n1_z < 0 || n1_z >= nz) return 0.0f;
        const float val_n1 = F[n1_x * stride_x + n1_y * stride_y + n1_z];
        if (std::abs(val_n1 - val) > epsilon) return 0.0f;
        const int n2_x = x + dx * 2, n2_y = y + dy * 2, n2_z = z + dz * 2;
        if (n2_x < 0 || n2_x >= nx || n2_y < 0 || n2_y >= ny || n2_z < 0 || n2_z >= nz) return (float)dx;
        const float val_n2 = F[n2_x * stride_x + n2_y * stride_y + n2_z];
        if (val_n2 < val - epsilon) return (float)dx;
        return 0.0f;
    };
    float sx = axis_shift(1, 0, 0);
    if (sx == 0.0f) sx = axis_shift(-1, 0, 0);
    float sy = axis_shift(0, 1, 0);
    if (sy == 0.0f) sy = axis_shift(0, -1, 0);
    float sz = axis_shift(0, 0, 1);
    if (sz == 0.0f) sz = axis_shift(0, 0, -1);
    return Eigen::Vector3f(sx, sy, sz);
}

/**
 * @brief Sub-voxel shift of a local maximum of a stored field (see the templated shift_voxel_center()).
 * @param field Voxel grid field.
 * @param idx Voxel of the maximum.
 * @param epsilon Tolerance for equal values.
 * @return Shift direction per axis (each -1, 0 or 1).
 */
inline Eigen::Vector3f shift_voxel_center(
    const VoxelGrid<float>& field,
    const Eigen::Vector3i& idx,
    float epsilon = 1e-5f)
{
    return shift_voxel_center(GridField{field.data.data()}, static_cast<int>(field.nx()), static_cast<int>(field.ny()),
                              static_cast<int>(field.nz()), idx, epsilon);
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

    size_t target_vol = 0;
    size_t recon_vol = 0;
    size_t intersection = 0;
    const size_t n = target.data.size();

    #pragma omp parallel for reduction(+:target_vol, recon_vol, intersection)
    for (size_t i = 0; i < n; ++i) {
        bool t = (target.data[i] > static_cast<T>(0));
        bool r = (reconstruction.data[i] > static_cast<U>(0));
        
        if (t) target_vol++;
        if (r) recon_vol++;
        if (t && r) intersection++;
    }

    if (target_vol == 0 && recon_vol == 0) return 1.0;
    if (target_vol == 0 || recon_vol == 0) return 0.0;
    
    return (2.0 * static_cast<double>(intersection)) / static_cast<double>(target_vol + recon_vol);
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
    T fill_value = static_cast<T>(1), float coefficient = 0.0f)
{
    if (sphere_table.rows() == 0) return;
    #pragma omp parallel for
    for (int i = 0; i < sphere_table.rows(); ++i) {
        float cx = sphere_table(i, 0);
        float cy = sphere_table(i, 1);
        float cz = sphere_table(i, 2);
        float radius_vox = sphere_table(i, 3);
        if (coefficient != 0.0f) {
            radius_vox = std::sqrt (radius_vox) * coefficient;
        }
        grid.sphere_kernel(cx, cy, cz, radius_vox, fill_value);
    }
}

/**
 * @brief Converts a voxel grid into a blocky mesh by generating faces for exposed voxels.
 * @tparam T Voxel grid data type.
 * @param grid Input voxel grid.
 * @param threshold Threshold for voxel occupancy.
 * @return STLMesh mesh structure.
 */
template <typename T>
STLMesh grid_to_mesh(const VoxelGrid<T>& grid, T threshold = static_cast<T>(0)) {
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

    STLMesh mesh;
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
