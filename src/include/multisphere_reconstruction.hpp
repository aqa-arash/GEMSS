#ifndef MULTISPHERE_RECONSTRUCTION_HPP
#define MULTISPHERE_RECONSTRUCTION_HPP

/**
 * @file multisphere_reconstruction.hpp
 * @brief Main reconstruction algorithms for multisphere-cpp.
 *
 * Provides multisphere reconstruction from voxel grids and triangle meshes.
 * Includes iterative peak detection, filtering, and conversion to physical units.
 *
 * @author Arash Moradian
 * @date 2026-03-09
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <type_traits>
#include <Eigen/Dense>
#ifdef HAVE_OPENMP
    #include <omp.h>
#endif

#include "multisphere_datatypes.hpp"
#include "multisphere_reconstruction_helpers.hpp"
#include "multisphere_io.hpp"
#include "multisphere_mesh_handler.hpp"
#include "multisphere_physics_computation.hpp"
#include "multisphere_config.hpp"

namespace MSS {

/**
 * @brief Construct multisphere from voxel grid.
 * Iteratively detects peaks, filters, and converts to physical units.
 * @tparam T VoxelGrid data type.
 * @param input_grid Input voxel grid.
 * @param config MultisphereConfig struct containing all configuration parameters:
 * @param min_center_distance_vox Minimum center distance (voxels).
 * @param min_radius_vox Minimum radius (voxels).
 * @param precision_target Target precision.
 * @param max_spheres Maximum number of spheres.
 * @param show_progress Show progress output.
 * @param sphere_table_in Optional initial sphere table.
 * @param compute_physics Whether to compute physical properties of the multisphere union.
 * @return SpherePack reconstruction result.
 */
template <typename T>
SpherePack multisphere_from_voxels(
    const VoxelGrid<T>& input_grid,
    const MultisphereConfig& config = MultisphereConfig()
) {


    
    VoxelGrid<float> original_distance(input_grid.nx(), input_grid.ny(), input_grid.nz(), input_grid.voxel_size, input_grid.origin);
    VoxelGrid<uint8_t> voxel_grid(input_grid.nx(), input_grid.ny(), input_grid.nz(), input_grid.voxel_size, input_grid.origin);
    Eigen::MatrixX4f sphere_table(0, 4);

    if (config.initial_sphere_table.has_value()) {
        sphere_table = *config.initial_sphere_table;
        std::cout << "Using provided initial sphere table with " << sphere_table.rows() << " spheres." << std::endl;
    } else {
        std::cout << "No initial sphere table provided. Starting with empty table." << std::endl;
    }

    
    if (config.min_center_distance_vox <= 1) {
        throw std::invalid_argument("min_center_distance_vox must be greater than 1 to avoid numerical issues.");
    }

    // Convert to uint8_t grid for distance transform if needed
    if constexpr (std::is_same<T, uint8_t>::value) {
        voxel_grid = input_grid;
    } else {
        #pragma omp parallel for
        for (size_t i = 0; i < input_grid.data.size(); ++i) {
            voxel_grid.data[i] = (input_grid.data[i] > static_cast<T>(0)) ? static_cast<T>(1) : static_cast<T>(0);
        }
    }

    original_distance = voxel_grid.distance_transform();

    #ifdef MULTISPHERE_DEBUG
    // [DEBUG] Distance field range
    float min_d = 1e9, max_d = -1e9;
    for (float v : original_distance.data) {
        if (v < min_d) min_d = v;
        if (v > max_d) max_d = v;
    }
    std::cout << "[DEBUG] Distance Field Range: [" << min_d << ", " << max_d << "]" << std::endl;
    if (max_d <= 0.0f) {
        std::cerr << "[ERROR] Distance field is zero/negative. Check EDT logic!" << std::endl;
    }
    #endif

    VoxelGrid<uint8_t> recon_mask(voxel_grid.nx(), voxel_grid.ny(), voxel_grid.nz(), voxel_grid.voxel_size, voxel_grid.origin);

    // iterative solver to compute sphere table
    Eigen::MatrixX4f final_sphere_table = compute_sphere_table(
        original_distance,
        voxel_grid,
        recon_mask,
        sphere_table,
        config.min_center_distance_vox,
        config.min_radius_vox,
        config.precision_target,
        config.max_spheres,
        config.show_progress
        );

    #ifdef MULTISPHERE_DEBUG
        // [DEBUG] Final spheres
        std::cout << "Final Spheres: " << final_sphere_table.rows() << std::endl;
        std::cout << " Sample Sphere Radii: ";
        for (int i = 0; i < std::min(5, int(final_sphere_table.rows())); ++i) {
            std::cout << "x = " << final_sphere_table(i, 0) << ", y = " << final_sphere_table(i, 1) << ", z = " << final_sphere_table(i, 2) << ", r = " << final_sphere_table(i, 3) << "\n";
        }
        std::cout << std::endl;
    #endif

    // =========================================================================
    // Prune Isolated Networks
    // =========================================================================
    if (config.prune_isolated_spheres && final_sphere_table.rows() > 1) {
        int original_count = final_sphere_table.rows();
        final_sphere_table = filter_largest_sphere_network(final_sphere_table);
        
        if (final_sphere_table.rows() < original_count) {
            if (config.show_progress) {
                std::cout << "Pruned " << (original_count - final_sphere_table.rows()) 
                          << " isolated spheres. Retained " << final_sphere_table.rows() 
                          << " in the primary network." << std::endl;
            }
            
            // Rebuild the mask so precision & physics calculations 
            // strictly reflect the pruned primary network.
            std::fill(recon_mask.data.begin(), recon_mask.data.end(), static_cast<uint8_t>(0));
            spheres_to_grid<uint8_t>(recon_mask, final_sphere_table);
        }
    }
    // =========================================================================

    float final_precision = compute_voxel_precision(voxel_grid, recon_mask);

    if(config.show_progress) {
        std::cout << "Final Precision: " << final_precision << std::endl;
    }

    // Convert Table to SpherePack (Physical units)
    Eigen::MatrixX3f centers_phys(final_sphere_table.rows(), 3);
    Eigen::VectorXf radii_phys(final_sphere_table.rows());

    for (int i = 0; i < final_sphere_table.rows(); ++i) {
        Eigen::Vector3f pos_vox = final_sphere_table.block<1, 3>(i, 0).transpose();
        centers_phys.row(i) = voxel_grid.origin + (pos_vox.array() + 0.5f).matrix() * voxel_grid.voxel_size;
        radii_phys(i) = (final_sphere_table(i, 3)) * voxel_grid.voxel_size;
    }

    #ifdef MULTISPHERE_DEBUG
        // [DEBUG] Final sphere sample (physical units)
        std::cout << "[DEBUG] Final Sphere Sample (Physical Units):" << std::endl;
        for (int i = 0; i < std::min(5, int(final_sphere_table.rows())); ++i) {
            std::cout << "x = " << centers_phys(i, 0) << ", y = " << centers_phys(i, 1) << ", z = " << centers_phys(i, 2) << ", r = " << radii_phys(i) << "\n";
        }
        std::cout << std::endl;
    #endif

    SpherePack result(centers_phys, radii_phys);
    result.precision = final_precision;
    if (config.compute_physics) {
        // Compute physical properties of the multisphere union
        compute_multisphere_physics(result, recon_mask);
    }
    

    return result;
}

/**
 * @brief Construct a multisphere representation directly from a triangle mesh.
 * Logically equivalent to the Python version, optimized for C++.
 * @param mesh Input FastMesh.
 * @param config MultisphereConfig struct containing all configuration parameters:
 * @param div Voxel grid division (resolution).
 * @param padding Grid padding.
 * @param min_center_distance_vox Minimum center distance (voxels).
 * @param min_radius_vox Minimum radius (voxels).
 * @param precision Target precision.
 * @param max_spheres Maximum number of spheres.
 * @param show_progress Show progress output.
 * @param confine_mesh Confine spheres to mesh boundary.
 * @param sphere_table Optional initial sphere table.
 * @param compute_physics Whether to compute physical properties of the multisphere union.
 * @return SpherePack reconstruction result.
 */
SpherePack multisphere_from_mesh(
    const FastMesh& mesh,
    const MultisphereConfig& config = MultisphereConfig()
) {
    if (mesh.is_empty()) {
        throw std::runtime_error("Cannot reconstruct from an empty mesh.");
    }

    // 1. Convert to VoxelGrid
    VoxelGrid<bool> voxel_grid = mesh_to_binary_grid(mesh, config.div, config.padding);
    #ifdef MULTISPHERE_DEBUG    
        std::cout << "Voxel grid created from mesh: " << voxel_grid.nx() << "x" << voxel_grid.ny() << "x" << voxel_grid.nz() << std::endl;
    #endif

    SpherePack sp;
    if (config.initial_sphere_table.has_value()) {
        // Convert Table to SpherePack (Physical units)
        Eigen::MatrixX3f centers_vox(config.initial_sphere_table.value().rows(), 3);
        Eigen::VectorXf radii_vox(config.initial_sphere_table.value().rows());

        for (int i = 0; i < config.initial_sphere_table.value().rows(); ++i) {
            Eigen::Vector3f pos_phys = config.initial_sphere_table.value().block<1, 3>(i, 0).transpose();
            centers_vox.row(i) = ((pos_phys - voxel_grid.origin).array() / voxel_grid.voxel_size) - 0.5f;
            radii_vox(i) = (config.initial_sphere_table.value().operator()(i, 3)) / voxel_grid.voxel_size;
        }
        Eigen::MatrixX4f sphere_table_vox(config.initial_sphere_table.value().rows(), 4);
        sphere_table_vox.block(0, 0, config.initial_sphere_table.value().rows(), 3) = centers_vox;
        sphere_table_vox.col(3) = radii_vox;

        MSS::MultisphereConfig config_vox = config; // Make a copy to modify
        config_vox.initial_sphere_table = sphere_table_vox; // Pass the voxel-space sphere table to the config for reconstruction

        #ifdef MULTISPHERE_DEBUG
        std::cout << "Converted initial sphere table to voxel units for reconstruction." << std::endl;
        #endif

        sp = multisphere_from_voxels(
            voxel_grid,
            config_vox
        );
    } else {
        #ifdef MULTISPHERE_DEBUG
            std::cout << "No initial sphere table provided. Starting with an empty table." << std::endl;
        #endif
        sp = multisphere_from_voxels(
            voxel_grid,
            config
        );
    }

    // 4. Boundary Adjustment (Optional)
    if (config.confine_mesh) {
        constrain_radii_to_sdf(sp, mesh);
    }

    return sp;
}


} // namespace MSS

#endif // MULTISPHERE_RECONSTRUCTION_HPP
