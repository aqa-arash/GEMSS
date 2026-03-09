#ifndef MULTISPHERE_MESH_HANDLER_HPP
#define MULTISPHERE_MESH_HANDLER_HPP

#include <vector>
#include <cmath>
#include <iostream>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <algorithm>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "multisphere_datatypes.hpp"
// 1. Save current warning state and ignore specific warning
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfree-nonheap-object"
#include <igl/signed_distance.h>
#include <igl/fast_winding_number.h>
#pragma GCC diagnostic pop

// --- Robust Voxelizer (Generalized Winding Number) ---
// Best for "Dirty" Meshes (Holes, Self-Intersections, Internal Faces)
inline VoxelGrid<bool> mesh_to_binary_grid(const FastMesh& mesh, int div, int padding = 2) {
    if (mesh.vertices.rows() == 0) throw std::invalid_argument("Mesh is empty.");

    // 1. Setup Grid & Bounds
    // ------------------------------------------------
    Eigen::Vector3f min_v = mesh.vertices.colwise().minCoeff().transpose();
    Eigen::Vector3f max_v = mesh.vertices.colwise().maxCoeff().transpose();
    Eigen::Vector3f extents = max_v - min_v;
    
    // Handle flat meshes safely
    float min_extent = extents.minCoeff();
    if (min_extent <= 1e-6) min_extent = (extents.maxCoeff() > 0 ? extents.maxCoeff() : 1.0f);

    float voxel_size = min_extent / static_cast<float>(div);
    
    int nx = std::ceil(extents.x() / voxel_size) + 2 * padding;
    int ny = std::ceil(extents.y() / voxel_size) + 2 * padding;
    int nz = std::ceil(extents.z() / voxel_size) + 2 * padding;
    
    Eigen::Vector3f origin = min_v.cast<float>() - (float)padding * voxel_size * Eigen::Vector3f::Ones();

    std::cout << "[Voxelizer] Grid: " << nx << "x" << ny << "x" << nz 
              << " | Method: Generalized Winding Number (Robust)" << std::endl;

    VoxelGrid<bool> grid(nx, ny, nz, voxel_size, origin);

    // 2. Prepare Query Points (Voxel Centers)
    // ------------------------------------------------
    // libigl expects a Matrix of points to test.
    Eigen::MatrixXf queries(nx * ny * nz, 3);
    
    // We use OpenMP to fill the query matrix efficiently
    #pragma omp parallel for collapse(3)
    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int z = 0; z < nz; ++z) {
                int idx = x * (ny * nz) + y * nz + z;
                queries(idx, 0) = static_cast<float>(origin.x() + (x + 0.5) * voxel_size);
                queries(idx, 1) = static_cast<float>(origin.y() + (y + 0.5) * voxel_size);
                queries(idx, 2) = static_cast<float>(origin.z() + (z + 0.5) * voxel_size);
            }
        }
    }

    // 3. Compute Winding Number
    // ------------------------------------------------
    // This is the "Magic" step.
    // Result > 0.5 is mathematically "Inside".
    // Result < 0.5 is "Outside".
    // Because it integrates over the whole surface, small holes don't break it.
    
    Eigen::VectorXf wn;
    
    // igl::fast_winding_number computes the approximation for all points at once.
    // It creates an internal hierarchy (Octree) automatically for speed.
    // Inputs: (Vertices, Triangles, QueryPoints, OutputVector)
    igl::fast_winding_number(mesh.vertices, mesh.triangles, queries, wn);

    // 4. Threshold & Fill Grid
    // ------------------------------------------------
    #pragma omp parallel for
    for (int i = 0; i < wn.size(); ++i) {
        // We use abs() because sometimes triangles are inverted (backwards).
        // The magnitude tells us if we are inside.
        if (std::abs(wn(i)) >= 0.5f) {
            grid.data[i] = true;
        } else {
            grid.data[i] = false;
        }
    }

    return grid;
}

/**
 * @brief Constrains sphere radii so they do not exceed the distance to the nearest surface.
 * * @param pack The SpherePack to modify.
 */
void constrain_radii_to_sdf(SpherePack& pack, const FastMesh& mesh) {
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

#endif