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

namespace MSS {

/**
 * @brief Construct multisphere from voxel grid.
 * Iteratively detects peaks, filters, and converts to physical units.
 * @tparam T VoxelGrid data type.
 * @param input_grid Input voxel grid.
 * @param min_center_distance_vox Minimum center distance (voxels).
 * @param min_radius_vox Minimum radius (voxels).
 * @param precision_target Target precision.
 * @param max_spheres Maximum number of spheres.
 * @param show_progress Show progress output.
 * @param sphere_table_in Optional initial sphere table.
 * @return SpherePack reconstruction result.
 */
template <typename T>
SpherePack multisphere_from_voxels(
    const VoxelGrid<T>& input_grid,
    int min_center_distance_vox = 2,
    std::optional<int> min_radius_vox = std::nullopt,
    std::optional<float> precision_target = std::nullopt,
    std::optional<int> max_spheres = std::nullopt,
    bool show_progress = true,
    std::optional<Eigen::MatrixX4f> sphere_table_in = std::nullopt
) {


    
    VoxelGrid<float> original_distance(input_grid.nx(), input_grid.ny(), input_grid.nz(), input_grid.voxel_size, input_grid.origin);
    VoxelGrid<uint8_t> voxel_grid(input_grid.nx(), input_grid.ny(), input_grid.nz(), input_grid.voxel_size, input_grid.origin);
    Eigen::MatrixX4f sphere_table(0, 4);

    if (sphere_table_in.has_value()) {
        sphere_table = *sphere_table_in;
        std::cout << "Using provided initial sphere table with " << sphere_table.rows() << " spheres." << std::endl;
    } else {
        std::cout << "No initial sphere table provided. Starting with empty table." << std::endl;
    }

    bool check_min_center_distance = true;
    if (min_center_distance_vox <= 1) {
        std::cout << "[WARNING] min_center_distance_vox <= 1." << std::endl;
        min_center_distance_vox = 2;
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

    int max_iter = 3; // [DEBUG] Max iteration for residual scaling
    int iter = 1;
    float weight_factor = 1.0f; // [DEBUG] Initial weight factor for residual
    VoxelGrid<uint8_t> recon_mask(voxel_grid.nx(), voxel_grid.ny(), voxel_grid.nz(), voxel_grid.voxel_size, voxel_grid.origin);
    VoxelGrid<float> residual(voxel_grid.nx(), voxel_grid.ny(), voxel_grid.nz(), voxel_grid.voxel_size, voxel_grid.origin);

    int prev_count = 0;
    bool peaks_found = true;
    while (iter <= max_iter) {

        #ifdef MULTISPHERE_DEBUG
             std::cout << "Iteration " << iter << ": Current spheres = " << sphere_table.rows() << std::endl;
        #endif
        // Termination: Max Spheres
        if (max_spheres.has_value() && sphere_table.rows() >= *max_spheres) break;

        VoxelGrid<float> summed_field(voxel_grid.nx(), voxel_grid.ny(), voxel_grid.nz());

        if (sphere_table.rows() > 0) {
            if (peaks_found) {
                iter = 1; // Reset iteration if we found peaks in the last round
                spheres_to_grid<uint8_t>(
                    recon_mask,
                    sphere_table.block(prev_count, 0, sphere_table.rows() - prev_count, 4)
                );

                // Termination: Precision
                if (precision_target.has_value()) {
                    float current_precision = compute_voxel_precision(voxel_grid, recon_mask);
                    if (show_progress) {
                        std::cout << "Weight factor = " << weight_factor << " #spheres = " << sphere_table.rows() << " Precision = " << current_precision << std::endl;
                    }
                    if (current_precision >= *precision_target) {
                        if (show_progress) {
                            std::cout << "Reached target precision of " << *precision_target << ". Terminating." << std::endl;
                        }
                        break;
                    }
                }

                VoxelGrid<float> spheres_distance = recon_mask.distance_transform();
                residual = residual_distance_field(original_distance, spheres_distance);
            }

            // summed_field = distance_field + residual
            #pragma omp parallel for
            for (size_t i = 0; i < summed_field.data.size(); ++i) {
                summed_field.data[i] = original_distance.data[i] + residual.data[i] * weight_factor;
            }
        } else {
            summed_field.data = original_distance.data;
        }

        // 2. Peak Detection
        Eigen::MatrixX4f peaks = peak_local_max_3d(summed_field, original_distance, min_center_distance_vox, min_radius_vox.value_or(1));

        // 3. Filter Peaks
        peaks = filter_and_shift_peaks(peaks, sphere_table, original_distance, min_center_distance_vox);
        if (peaks.rows() == 0) {
            #ifdef MULTISPHERE_DEBUG
                std::cout << " No peaks found after filtering, for iter = " << iter << std::endl;
            #endif
            if (iter >= max_iter) {
                std::cout << "No peaks found after filtering, terminating." << std::endl;
                break;
            } else {
                peaks_found = false;
                ++iter;
                ++weight_factor; // [DEBUG] Increase weight factor to encourage new peaks
                #ifdef MULTISPHERE_DEBUG
                    std::cout << " Continuing to next iteration." << std::endl;
                #endif
                continue;
            }
        }
        peaks_found = true;

        #ifdef MULTISPHERE_DEBUG
        // [DEBUG] Peaks after filtering
        std::cout << " Peaks after Filtering: " << peaks.rows() << std::endl;
        std::cout << " Sample Peak Radii after Filtering: ";
        for (int i = 0; i < std::min(5, int(peaks.rows())); ++i) {
            std::cout << "x = " << peaks(i, 0) << ", y = " << peaks(i, 1) << ", z = " << peaks(i, 2) << ", r = " << peaks(i, 3) << "\n";
        }
        std::cout << std::endl;
        #endif
        // 4. Append Spheres
        prev_count = sphere_table.rows();
        append_sphere_table(sphere_table, peaks, max_spheres);

        if (min_radius_vox.has_value() && sphere_table.rows() == prev_count) {
            #ifdef MULTISPHERE_DEBUG
                std::cout << " No new spheres added that meet min_radius_vox. Terminating." << std::endl;
            #endif
            break;
        }
    }
    #ifdef MULTISPHERE_DEBUG
        // [DEBUG] Final spheres
        std::cout << "Final Spheres: " << sphere_table.rows() << std::endl;
        std::cout << " Sample Sphere Radii: ";
        for (int i = 0; i < std::min(5, int(sphere_table.rows())); ++i) {
            std::cout << "x = " << sphere_table(i, 0) << ", y = " << sphere_table(i, 1) << ", z = " << sphere_table(i, 2) << ", r = " << sphere_table(i, 3) << "\n";
        }
        std::cout << std::endl;
    #endif
    // Convert Table to SpherePack (Physical units)
    Eigen::MatrixX3f centers_phys(sphere_table.rows(), 3);
    Eigen::VectorXf radii_phys(sphere_table.rows());

    for (int i = 0; i < sphere_table.rows(); ++i) {
        Eigen::Vector3f pos_vox = sphere_table.block<1, 3>(i, 0).transpose();
        centers_phys.row(i) = voxel_grid.origin + (pos_vox.array() + 0.5f).matrix() * voxel_grid.voxel_size;
        radii_phys(i) = (sphere_table(i, 3)) * voxel_grid.voxel_size;
    }

    #ifdef MULTISPHERE_DEBUG
        // [DEBUG] Final sphere sample (physical units)
        std::cout << "[DEBUG] Final Sphere Sample (Physical Units):" << std::endl;
        for (int i = 0; i < std::min(5, int(sphere_table.rows())); ++i) {
            std::cout << "x = " << centers_phys(i, 0) << ", y = " << centers_phys(i, 1) << ", z = " << centers_phys(i, 2) << ", r = " << radii_phys(i) << "\n";
        }
        std::cout << std::endl;
    #endif
    return SpherePack(centers_phys, radii_phys);
}

/**
 * @brief Construct a multisphere representation directly from a triangle mesh.
 * Logically equivalent to the Python version, optimized for C++.
 * @param mesh Input FastMesh.
 * @param div Voxel grid division (resolution).
 * @param padding Grid padding.
 * @param min_center_distance_vox Minimum center distance (voxels).
 * @param min_radius_vox Minimum radius (voxels).
 * @param precision Target precision.
 * @param max_spheres Maximum number of spheres.
 * @param show_progress Show progress output.
 * @param confine_mesh Confine spheres to mesh boundary.
 * @param sphere_table Optional initial sphere table.
 * @return SpherePack reconstruction result.
 */
SpherePack multisphere_from_mesh(
    const FastMesh& mesh,
    int div = 100,
    int padding = 2,
    int min_center_distance_vox = 4,
    std::optional<int> min_radius_vox = std::nullopt,
    std::optional<float> precision = std::nullopt,
    std::optional<int> max_spheres = std::nullopt,
    bool show_progress = true,
    bool confine_mesh = false,
    std::optional<Eigen::MatrixX4f> sphere_table = std::nullopt
) {
    if (mesh.is_empty()) {
        throw std::runtime_error("Cannot reconstruct from an empty mesh.");
    }

    // 1. Convert to VoxelGrid
    VoxelGrid<bool> voxel_grid = mesh_to_binary_grid(mesh, div, padding);

    SpherePack sp;
    if (sphere_table.has_value()) {
        // Convert Table to SpherePack (Physical units)
        Eigen::MatrixX3f centers_vox(sphere_table->rows(), 3);
        Eigen::VectorXf radii_vox(sphere_table->rows());

        for (int i = 0; i < sphere_table->rows(); ++i) {
            Eigen::Vector3f pos_phys = sphere_table->block<1, 3>(i, 0).transpose();
            centers_vox.row(i) = ((pos_phys - voxel_grid.origin).array() / voxel_grid.voxel_size) - 0.5f;
            radii_vox(i) = (sphere_table->operator()(i, 3)) / voxel_grid.voxel_size;
        }
        Eigen::MatrixX4f sphere_table_vox(sphere_table->rows(), 4);
        sphere_table_vox.block(0, 0, sphere_table->rows(), 3) = centers_vox;
        sphere_table_vox.col(3) = radii_vox;

        sp = multisphere_from_voxels(
            voxel_grid,
            min_center_distance_vox,
            min_radius_vox,
            precision,
            max_spheres,
            show_progress,
            sphere_table_vox
        );
    } else {
        #ifdef MULTISPHERE_DEBUG
            std::cout << "No initial sphere table provided. Starting with an empty table." << std::endl;
        #endif
        sp = multisphere_from_voxels(
            voxel_grid,
            min_center_distance_vox,
            min_radius_vox,
            precision,
            max_spheres,
            show_progress
        );
    }

    // 4. Boundary Adjustment (Optional)
    if (confine_mesh) {
        constrain_radii_to_sdf(sp, mesh);
    }

    return sp;
}

#endif // MULTISPHERE_RECONSTRUCTION_HPP

} // namespace MSS