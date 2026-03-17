#ifndef MULTISPHERE_MESH_HANDLER_HPP
#define MULTISPHERE_MESH_HANDLER_HPP

/**
 * @file multisphere_mesh_handler.hpp
 * @brief Mesh voxelization and sphere pack constraint utilities for multisphere-cpp.
 *
 * Provides robust mesh-to-voxel conversion and sphere radius constraint based on signed distance.
 *
 * @author Arash Moradian
 * @date 2026-03-09
 */

#include <vector>
#include <cmath>
#include <iostream>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <algorithm>

#ifdef MULTISPHERE_DEBUG
    #include <chrono>
#endif

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "multisphere_datatypes.hpp"

// Suppress specific warnings for libigl
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfree-nonheap-object"
#include <igl/signed_distance.h>
#include <igl/fast_winding_number.h>
#pragma GCC diagnostic pop

// Namespace MSS for all library code
namespace MSS {

/** 
 * @brief Get the minimum axis-aligned bounding box dimension.
 * @param mesh Input mesh.
 * @return Minimum extent of the AABB.
 * Note: This is used to determine voxel size for grid construction. Where h = min_extent / div.
 */

float get_min_AABB(const FastMesh & mesh)  {
    // 1. Setup Grid & Bounds
    Eigen::Vector3f min_v = mesh.vertices.colwise().minCoeff().transpose();
    Eigen::Vector3f max_v = mesh.vertices.colwise().maxCoeff().transpose();
    Eigen::Vector3f extents = max_v - min_v;

    float min_extent = extents.minCoeff();
    if (min_extent <= 1e-6) min_extent = (extents.maxCoeff() > 0 ? extents.maxCoeff() : 1.0f);

    return min_extent;
}

/**
 * @brief Robust voxelizer using generalized winding number.
 *        Best for "dirty" meshes (holes, self-intersections, internal faces).
 * @param mesh Input FastMesh.
 * @param div Voxel grid division (resolution).
 * @param padding Grid padding.
 * @return VoxelGrid<uint8_t> representing mesh occupancy.
 */
inline VoxelGrid<uint8_t> mesh_to_binary_grid(const FastMesh& mesh, int div, int padding = 2) {
    if (mesh.vertices.rows() == 0) throw std::invalid_argument("Mesh is empty.");
    // 1. Setup Grid & Bounds
    Eigen::Vector3f min_v = mesh.vertices.colwise().minCoeff().transpose();
    Eigen::Vector3f max_v = mesh.vertices.colwise().maxCoeff().transpose();
    Eigen::Vector3f extents = max_v - min_v;
    float min_extent = extents.minCoeff();
    if (min_extent <= 1e-6) min_extent = (extents.maxCoeff() > 0 ? extents.maxCoeff() : 1.0f);
    float voxel_size = min_extent / static_cast<float>(div);

    int nx = std::ceil(extents.x() / voxel_size) + 2 * padding;
    int ny = std::ceil(extents.y() / voxel_size) + 2 * padding;
    int nz = std::ceil(extents.z() / voxel_size) + 2 * padding;

    Eigen::Vector3f origin = min_v.cast<float>() - static_cast<float>(padding) * voxel_size * Eigen::Vector3f::Ones();
    #ifdef MULTISPHERE_DEBUG
        std::cout << "[Voxelizer] Grid: " << nx << "x" << ny << "x" << nz 
                << " | Method: Generalized Winding Number (Robust)" << std::endl;
    #endif
    VoxelGrid<uint8_t> grid(nx, ny, nz, voxel_size, origin);

    // 2. Prepare Query Points (Voxel Centers)
    Eigen::MatrixXf queries(nx * ny * nz, 3);
    #ifdef MULTISPHERE_DEBUG
        auto start_time = std::chrono::high_resolution_clock::now();
    #endif
    #pragma omp parallel for collapse(2)
    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int z = 0; z < nz; ++z) {
                int idx = x * (ny * nz) + y * nz + z;
                queries(idx, 0) = static_cast<float>(origin.x() + (x + 0.5f) * voxel_size);
                queries(idx, 1) = static_cast<float>(origin.y() + (y + 0.5f) * voxel_size);
                queries(idx, 2) = static_cast<float>(origin.z() + (z + 0.5f) * voxel_size);
            }
        }
    }

    #ifdef MULTISPHERE_DEBUG
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        std::cout << "[Voxelizer] Query points generated in " << elapsed.count() << " seconds." << std::endl;
    #endif
    
    #ifdef MULTISPHERE_DEBUG
        start_time = std::chrono::high_resolution_clock::now();
    #endif
    // 3. Compute Winding Number
    Eigen::VectorXf wn;
    igl::fast_winding_number(mesh.vertices, mesh.triangles, queries, wn);

    #ifdef MULTISPHERE_DEBUG
        end_time = std::chrono::high_resolution_clock::now();
        elapsed = end_time - start_time;
        std::cout << "[Voxelizer] Winding number computed in " << elapsed.count() << " seconds." << std::endl;
    #endif

    #ifdef MULTISPHERE_DEBUG
        std::cout << "[Voxelizer] Winding number computed for " << queries.rows() << " query points." << std::endl;
    #endif
    // 4. Threshold & Fill Grid
    #pragma omp parallel for
    for (int i = 0; i < wn.size(); ++i) {
        // Use abs() for robustness against inverted triangles
        grid.data[i] = (std::abs(wn(i)) >= 0.5f);
    }

    return grid;
}

/**
 * @brief Constrains sphere radii so they do not exceed the distance to the nearest surface.
 * @param pack SpherePack to modify.
 * @param mesh FastMesh reference mesh.
 */
inline void constrain_radii_to_sdf(SpherePack& pack, const FastMesh& mesh) {
    if (mesh.is_empty()) return;

    Eigen::VectorXf sdf;
    Eigen::VectorXi I;
    Eigen::MatrixXf C, N;

    // Compute signed distance
    igl::signed_distance(
        pack.centers, 
        mesh.vertices,
        mesh.triangles, 
        igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, 
        sdf, I, C, N
    );

    // Update radii: snap to surface (abs(sdf))
    pack.radii = pack.radii.array().min(sdf.array().abs());
}

} // namespace MSS

#endif // MULTISPHERE_MESH_HANDLER_HPP