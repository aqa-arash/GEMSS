#ifndef GEMSS_RECONSTRUCTION_HPP
#define GEMSS_RECONSTRUCTION_HPP


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
#include <Eigen/Core>
#ifdef HAVE_OPENMP
    #include <omp.h>
#endif

#include "GEMSS_datatypes.hpp"
#include "GEMSS_reconstruction_helpers.hpp"
#include "GEMSS_io.hpp"
#include "GEMSS_mesh_handler.hpp"
#include "GEMSS_physics_computation.hpp"
#include "GEMSS_config.hpp"

namespace GEMSS {

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
inline SpherePack multisphere_from_voxels(
    const VoxelGrid<uint8_t>& input_grid,
    const MultisphereConfig& config = MultisphereConfig()
) {


    
    VoxelGrid<float> original_distance(input_grid.nx(), input_grid.ny(), input_grid.nz(), input_grid.voxel_size, input_grid.origin);
    Eigen::MatrixX4f sphere_table(0, 4);

    
    if (config.search_window <= 1) {
        throw std::invalid_argument("search_window must be greater than 1 to avoid numerical issues.");
    }

 
    original_distance = input_grid.distance_transform();

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

    VoxelGrid<uint8_t> recon_mask(input_grid.nx(), input_grid.ny(), input_grid.nz(), input_grid.voxel_size, input_grid.origin);

    // iterative solver to compute sphere table
    Eigen::MatrixX4f final_sphere_table = compute_sphere_table(
        original_distance,
        input_grid,
        recon_mask,
        config
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

    float final_precision = compute_voxel_precision(input_grid, recon_mask);

    if(config.show_progress) {
        std::cout << "Final Precision: " << final_precision << std::endl;
    }

    // Convert Table to SpherePack (Physical units)
    Eigen::MatrixX3f centers_phys(final_sphere_table.rows(), 3);
    Eigen::VectorXf radii_phys(final_sphere_table.rows());

    for (int i = 0; i < final_sphere_table.rows(); ++i) {
        Eigen::Vector3f pos_vox = final_sphere_table.block<1, 3>(i, 0).transpose();
        centers_phys.row(i) = input_grid.origin + (pos_vox.array() + 0.5f).matrix() * input_grid.voxel_size;
        radii_phys(i) = (final_sphere_table(i, 3)) * input_grid.voxel_size;
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
    if (config.compute_physics == 1) {
        // Compute physical properties of the multisphere union
        compute_multisphere_physics(result, recon_mask);
    }
    if (config.compute_physics == 2) {
        // Compute physical properties based on original mesh (if available)
        compute_multisphere_physics(result, input_grid);
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
inline SpherePack multisphere_from_mesh(
    const FastMesh& mesh,
    const MultisphereConfig& config = MultisphereConfig()
) {
    if (mesh.is_empty()) {
        throw std::runtime_error("Cannot reconstruct from an empty mesh.");
    }

    GEMSS::MultisphereConfig config_vox = config; // Make a copy to modify


    // 1. Convert to VoxelGrid
    VoxelGrid<uint8_t> voxel_grid = mesh_to_binary_grid(mesh, config_vox);


    #ifdef MULTISPHERE_DEBUG    
        std::cout << "Voxel grid created from mesh: " << voxel_grid.nx() << "x" << voxel_grid.ny() << "x" << voxel_grid.nz() << std::endl;
    #endif

    // shift and scale initial sphere table from physical units to voxel units if provided
    if (config.initial_sphere_table.rows() > 0) {
        // Copy the entire table to start
        config_vox.initial_sphere_table = config.initial_sphere_table;

        // 1. Shift centers by origin (broadcasting the origin row-wise)
        config_vox.initial_sphere_table.leftCols<3>().rowwise() -= voxel_grid.origin.transpose();

        // 2. Scale centers by voxel size and offset by 0.5
        config_vox.initial_sphere_table.leftCols<3>().array() = 
            (config_vox.initial_sphere_table.leftCols<3>().array() / voxel_grid.voxel_size) - 0.5f;

        // 3. Scale radii by voxel size
        config_vox.initial_sphere_table.col(3).array() /= voxel_grid.voxel_size;

        #ifdef MULTISPHERE_DEBUG
        std::cout << "Converted initial sphere table to voxel units for reconstruction." << std::endl;
        #endif
    }

    SpherePack sp;
    sp = multisphere_from_voxels(
            voxel_grid,
            config_vox
        );

    // 4. Boundary Adjustment (Optional)
    if (config.confine_mesh) {
        constrain_radii_to_sdf(sp, mesh);
    }

    return sp;
}


} // namespace MSS
#endif // GEMSS_RECONSTRUCTION_HPP
