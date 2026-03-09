#ifndef MULTISPHERE_VOXEL_PROCESSING_HPP
#define MULTISPHERE_VOXEL_PROCESSING_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <Eigen/Dense>
#ifdef HAVE_OPENMP
    #include <omp.h>
#endif

#include "multisphere_datatypes.hpp"

// --- Apply Kernel to Voxel Grid ---
enum class KernelMode { ADD, SUBTRACT, ZERO };

template<typename T, typename K>
void apply_kernel_to_grid(
    VoxelGrid<T>& grid,
    const Eigen::Vector3f& center,
    const VoxelGrid<K>& kernel, // Using VoxelGrid as a container for kernel data
    KernelMode mode = KernelMode::ADD,
    float scale = 1.0f) 
{
    // 1. Boundary & Dimension Checks
    const int dim = 3; 
    Eigen::Vector3i k_shape(kernel.nx(), kernel.ny(), kernel.nz());
    Eigen::Vector3i v_shape(grid.nx(), grid.ny(), grid.nz());

    // Parity check (Even kernel -> .5 coords, Odd kernel -> integer coords)
    for (int i = 0; i < dim; ++i) {
        double floor_val;
        double frac = std::modf(center[i], &floor_val);
        if (k_shape[i] % 2 == 0) {
            if (std::abs(frac - 0.5) > 1e-7) throw std::invalid_argument("Even kernel requires .5 center.");
        } else {
            if (std::abs(frac) > 1e-7) throw std::invalid_argument("Odd kernel requires integer center.");
        }
    }

    // 2. Calculate Overlap Region (llf: lower-left-front, trr: top-right-rear)
    Eigen::Vector3f half = (k_shape.cast<float>().array() - 1.0f) / 2.0f;
    Eigen::Vector3i llf_k = (center - half).array().round().cast<int>();
    Eigen::Vector3i trr_k = (center + half).array().round().cast<int>();

    // Clipping against volume boundaries
    Eigen::Vector3i start_v = llf_k.cwiseMax(0);
    Eigen::Vector3i end_v = trr_k.cwiseMin((v_shape.array() - 1).matrix()); // Add .matrix()
    
    // Calculate kernel internal offsets if clipped
    Eigen::Vector3i k_offset = start_v - llf_k;

    // 3. Execution (Cache-friendly nested loops)
    // We use OpenMP for parallelization if the kernel is large
    #pragma omp parallel for collapse(2) if((end_v[0] - start_v[0]) > 32)
    for (int x = start_v[0]; x <= end_v[0]; ++x) {
        for (int y = start_v[1]; y <= end_v[1]; ++y) {
            // Inner loop is contiguous in memory (Z-axis) for SIMD optimization
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


// --- Helper: Check One Axis Direction ---
// Returns the shift value (e.g. 1.0 or -1.0) if a 2-voxel peak is detected.
// Returns 0.0 if it is a standard integer peak or a wide plateau.
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

    // 1. Check Immediate Neighbor (+1 step)
    int n1_x = cx + dx;
    int n1_y = cy + dy;
    int n1_z = cz + dz;

    // Bounds check
    if (n1_x < 0 || n1_x >= nx || n1_y < 0 || n1_y >= ny || n1_z < 0 || n1_z >= nz) return 0.0f;

    float val_n1 = field(n1_x, n1_y, n1_z);

    // If neighbor is not roughly equal, it's not a candidate for shifting.
    if (std::abs(val_n1 - current_val) > epsilon) return 0.0f;

    // 2. Check "Neighbor of Neighbor" (+2 steps) -> The Plateau Check
    int n2_x = cx + (dx * 2);
    int n2_y = cy + (dy * 2);
    int n2_z = cz + (dz * 2);

    // If we hit a wall immediately after the equal neighbor, treat it as a valid sub-voxel peak.
    if (n2_x < 0 || n2_x >= nx || n2_y < 0 || n2_y >= ny || n2_z < 0 || n2_z >= nz) return (float)dx;

    float val_n2 = field(n2_x, n2_y, n2_z);

    // [DECISION LOGIC]
    // A. If the value drops off after 1 neighbor, the peak is exactly 2 voxels wide.
    //    This means the true center is between current and n1. -> SHIFT.
    if (val_n2 < current_val - epsilon) {
        return (float)dx;
    }

    // B. If the value stays equal (or rises), the feature is >2 voxels wide (Plateau).
    //    Shifting 0.5 is meaningless here. -> NO SHIFT.
    return 0.0f; 
}


/**
 * Checks neighbors to see if a local maximum is centered on a voxel 
 * or between them (half-voxel shift).
 */
// --- Main Shift Function ---
inline Eigen::Vector3f shift_voxel_center(
    const VoxelGrid<float>& field, 
    const Eigen::Vector3i& idx, 
    float epsilon=1e-5) 
{
    int x = idx.x(), y = idx.y(), z = idx.z();
    float val = field(x, y, z);

    // Check X (+1 then -1)
    float sx = check_axis_shift(field, x, y, z, 1, 0, 0, val, epsilon);
    if (sx == 0.0f) sx = check_axis_shift(field, x, y, z, -1, 0, 0, val, epsilon);

    // Check Y (+1 then -1)
    float sy = check_axis_shift(field, x, y, z, 0, 1, 0, val, epsilon);
    if (sy == 0.0f) sy = check_axis_shift(field, x, y, z, 0, -1, 0, val, epsilon);

    // Check Z (+1 then -1)
    float sz = check_axis_shift(field, x, y, z, 0, 0, 1, val, epsilon);
    if (sz == 0.0f) sz = check_axis_shift(field, x, y, z, 0, 0, -1, val, epsilon);

    return Eigen::Vector3f(sx, sy, sz);
}


/**
 * Computes precision based on mismatch between target and reconstruction.
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

    // OMP reduction allows threads to sum their local counts then merge them
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


//Note: this might be wrong 
//TODO: double check the logic 
/** 
 * Renders spheres into a VoxelGrid. 
 * sphere_table columns: [x, y, z, diameter_vox]
 */
template <typename T>
void spheres_to_grid( VoxelGrid<T>& grid,
    const Eigen::MatrixX4f& sphere_table,
    T fill_value = static_cast<T>(1)) 
{
    if (sphere_table.rows() == 0) return;

    // We loop through spheres. For each sphere, apply_kernel_to_grid is called.
    // Note: If spheres are large, apply_kernel_to_grid internally uses OMP.
    #pragma omp parallel for
    for (int i = 0; i < sphere_table.rows(); ++i) {
        float cx = sphere_table(i, 0);
        float cy = sphere_table(i, 1);
        float cz = sphere_table(i, 2);
        float radius_vox = sphere_table(i, 3);
        
        if (radius_vox <= 0) continue;
        grid.sphere_kernel(cx, cy, cz, radius_vox, fill_value);


    }
}



/**
 * Converts a VoxelGrid into a "blocky" mesh by generating faces for exposed voxels.
 * This result can be passed directly to save_mesh_to_stl().
 */
template <typename T>
FastMesh grid_to_mesh(const VoxelGrid<T>& grid, T threshold = static_cast<T>(0)) {
    std::vector<Eigen::Vector3f> out_verts;
    std::vector<Eigen::Vector3i> out_tris;

    // Estimated reservation: Assuming roughly 10% surface area to volume ratio
    out_verts.reserve(grid.data.size() / 5);
    out_tris.reserve(grid.data.size() / 5);

    int nx = grid.nx();
    int ny = grid.ny();
    int nz = grid.nz();
    float vs = grid.voxel_size;

    // Direction offsets: +x, -x, +y, -y, +z, -z
    const int dx[6] = {1, -1, 0, 0, 0, 0};
    const int dy[6] = {0, 0, 1, -1, 0, 0};
    const int dz[6] = {0, 0, 0, 0, 1, -1};

    // Vertex offsets for the 6 faces (Quad corners relative to 0,0,0)
    // Order ensures Counter-Clockwise (CCW) winding for normals pointing outward
    const float face_verts[6][4][3] = {
        {{1,0,0}, {1,1,0}, {1,1,1}, {1,0,1}}, // +x (Right)
        {{0,0,0}, {0,0,1}, {0,1,1}, {0,1,0}}, // -x (Left)
        {{0,1,0}, {0,1,1}, {1,1,1}, {1,1,0}}, // +y (Top)
        {{0,0,0}, {1,0,0}, {1,0,1}, {0,0,1}}, // -y (Bottom)
        {{0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}}, // +z (Front)
        {{0,0,0}, {0,1,0}, {1,1,0}, {1,0,0}}  // -z (Back)
    };

    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int z = 0; z < nz; ++z) {
                // Skip empty voxels
                if (grid(x, y, z) <= threshold) continue;

                // Check all 6 neighbors
                for (int f = 0; f < 6; ++f) {
                    int nx_idx = x + dx[f];
                    int ny_idx = y + dy[f];
                    int nz_idx = z + dz[f];

                    bool draw_face = false;

                    // If neighbor is out of bounds, it's an edge -> Draw Face
                    if (nx_idx < 0 || nx_idx >= nx || 
                        ny_idx < 0 || ny_idx >= ny || 
                        nz_idx < 0 || nz_idx >= nz) {
                        draw_face = true;
                    } 
                    // If neighbor is empty, it's a surface -> Draw Face
                    else if (grid(nx_idx, ny_idx, nz_idx) <= threshold) {
                        draw_face = true;
                    }

                    if (draw_face) {
                        int base_idx = static_cast<int>(out_verts.size());

                        // Generate 4 vertices for this face
                        for (int v = 0; v < 4; ++v) {
                            Eigen::Vector3f pos;
                            pos.x() = grid.origin.x() + (x + face_verts[f][v][0]) * vs;
                            pos.y() = grid.origin.y() + (y + face_verts[f][v][1]) * vs;
                            pos.z() = grid.origin.z() + (z + face_verts[f][v][2]) * vs;
                            out_verts.push_back(pos);
                        }

                        // Add 2 Triangles (0-1-2, 0-2-3)
                        out_tris.push_back({base_idx, base_idx + 1, base_idx + 2});
                        out_tris.push_back({base_idx, base_idx + 2, base_idx + 3});
                    }
                }
            }
        }
    }

    // Convert to FastMesh format
    FastMesh mesh;
    mesh.vertices.resize(out_verts.size(), 3);
    mesh.triangles.resize(out_tris.size(), 3);

    // Using OMP for copy since vectors might be large
    #pragma omp parallel for
    for (size_t i = 0; i < out_verts.size(); ++i) mesh.vertices.row(i) = out_verts[i];
    #pragma omp parallel for
    for (size_t i = 0; i < out_tris.size(); ++i) mesh.triangles.row(i) = out_tris[i];

    return mesh;
}


#endif

