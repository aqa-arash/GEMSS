#ifndef MULTISPHERE_RECONSTRUCTION_HPP
#define MULTISPHERE_RECONSTRUCTION_HPP

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
#include "multisphere_utils.hpp"
#include "multisphere_io.hpp"
#include "multisphere_mesh_handler.hpp"


/**
 * Main Algorithm: Construct multisphere from voxel grid.
 */
template <typename T>
SpherePack multisphere_from_voxels(
    const VoxelGrid<T>& input_grid,
    int min_center_distance_vox = 2,
    int max_iter = 10,
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
    }
    else {
        std::cout << "No initial sphere table provided. Starting with empty table." << std::endl;
    }
    
    bool check_min_center_distance = true;
    if (min_center_distance_vox <= 1) {
        std::cout<<"[WARNING] min_center_distance_vox <= 1."<<std::endl;
        min_center_distance_vox = 2;
    }
    // check if T is uint8_t or bool for optimization
    if constexpr (std::is_same<T, uint8_t>::value) {
        voxel_grid = input_grid;
    } else {
        // Convert to uint8_t grid for distance transform
        #pragma omp parallel for
        for (size_t i = 0; i < input_grid.data.size(); ++i) {
            voxel_grid.data[i] = (input_grid.data[i] > static_cast<T>(0)) ? static_cast<T>(1) : static_cast<T>(0);
        }
    }
    // check type T and compute distance field if needed
    original_distance = voxel_grid.distance_transform();
    


    // [DEBUG] to verify field range
    float min_d = 1e9, max_d = -1e9;
    for(float v : original_distance.data) {
        if(v < min_d) min_d = v;
        if(v > max_d) max_d = v;
    }
    std::cout << "[DEBUG] Distance Field Range: [" << min_d << ", " << max_d << "]" << std::endl;
    if (max_d <= 0.0f) {
        std::cerr << "[ERROR] Distance field is zero/negative. Check EDT logic!" << std::endl;
    }
    // [END DEBUG]

    int iter = 1;

    VoxelGrid<uint8_t> recon_mask (voxel_grid.nx(),voxel_grid.ny(),voxel_grid.nz(), voxel_grid.voxel_size, voxel_grid.origin);
    
    VoxelGrid<float> residual(voxel_grid.nx(), voxel_grid.ny(), voxel_grid.nz(), voxel_grid.voxel_size, voxel_grid.origin);

    int prev_count = 0;
    bool peaks_found = true;
    while (true) {
        std::cout << "Iteration " << iter << ": Current spheres = " << sphere_table.rows() << std::endl;
        // Termination: Max Spheres
        if (max_spheres.has_value() && sphere_table.rows() >= *max_spheres) break;

        VoxelGrid<float> summed_field(voxel_grid.nx(), voxel_grid.ny(), voxel_grid.nz());

        if (sphere_table.rows() > 0) {

            if (peaks_found)  { 

                // uint8_t is thread-safe and acts as a boolean (0 vs 1).
                spheres_to_grid<uint8_t>( recon_mask, sphere_table.block(prev_count, 0, sphere_table.rows() - prev_count, 4) // new spheres only
                );

                // Termination: Precision
                if (precision_target.has_value()) {
                    float current_precision = compute_voxel_precision(voxel_grid, recon_mask);
                    std::cout<< " Iteration " << iter-1 << ": Current Precision = " << current_precision << std::endl;
                    if (current_precision >= *precision_target) {
                        std::cout << " Termination: Reached target precision of " << *precision_target << std::endl;
                        break;}
                    }

        
            VoxelGrid<float> spheres_distance = recon_mask.distance_transform();
            residual = residual_distance_field(original_distance, spheres_distance);
            }
            
                // summed_field = distance_field + residual
                #pragma omp parallel for
                for (size_t i = 0; i < summed_field.data.size(); ++i) {
                    summed_field.data[i] = original_distance.data[i] + residual.data[i] * iter ; 
                }
        } else {
            summed_field.data = original_distance.data;
        }

        // 2. Peak Detection
        Eigen::MatrixX4f peaks = peak_local_max_3d(summed_field, original_distance, min_center_distance_vox, min_radius_vox.value_or(1));
        
        /* [DEBUG]
        std::cout << " Detected Peaks: " << peaks.rows() << std::endl;
        std::cout<< " Sample Peak Radii: ";
        for (int i = 0; i < std::min(5, int(peaks.rows())); ++i) {
            std::cout <<"x = " << peaks(i, 0) << ", y = " << peaks(i, 1) << ", z = " << peaks(i, 2) << ", r = " << peaks(i, 3) << "\n";
        }
        std::cout << std::endl;
        */

        // 3. Filter Peaks
        peaks = filter_and_shift_peaks(peaks, sphere_table,original_distance, min_center_distance_vox);
        if (peaks.rows() == 0) {
            std::cout << " No peaks found after filtering, for iter = " << iter << std::endl;
            if (iter >= max_iter) {
                std::cout << " Terminating." << std::endl;
                break;
              }
            else {
                peaks_found = false;
                ++iter;
                std::cout << " Continuing to next iteration." << std::endl;
                continue;
            }
        }
        peaks_found = true;

        // [DEBUG]
        std::cout << " Peaks after Filtering: " << peaks.rows() << std::endl;
        std::cout<< " Sample Peak Radii after Filtering: ";
        for (int i = 0; i < std::min(5, int(peaks.rows())); ++i) {
            std::cout <<"x = " << peaks(i, 0) << ", y = " << peaks(i, 1) << ", z = " << peaks(i, 2) << ", r = " << peaks(i, 3) << "\n";
        }
        std::cout << std::endl;
        

        // 4. Append Spheres
        prev_count = sphere_table.rows();
        append_sphere_table(sphere_table, peaks, max_spheres);

        if (min_radius_vox.has_value() && sphere_table.rows() == prev_count) {
            std::cout << " No new spheres added that meet min_radius_vox. Terminating." << std::endl;
            break;}

    }


        // debug 
    std::cout << "Final Spheres: " << sphere_table.rows() << std::endl;
    std::cout<< " Sample Sphere Radii: ";
    for (int i = 0; i < std::min(5, int(sphere_table.rows())); ++i) {
        std::cout <<"x = " << sphere_table(i, 0) << ", y = " << sphere_table(i, 1) << ", z = " << sphere_table(i, 2) << ", r = " << sphere_table(i, 3) << "\n";
    }
    std::cout << std::endl;
 /* 
    //sort spheres by radius descending
    std::vector<std::pair<float, int>> radius_index_pairs;
    for (int i = 0; i < sphere_table.rows(); ++i) {
        radius_index_pairs.emplace_back(sphere_table(i, 3), i);
    }
    std::sort(radius_index_pairs.begin(), radius_index_pairs.end(), std::greater<std::pair<float, int>>());
    Eigen::MatrixX4f sorted_sphere_table(sphere_table.rows(), 4);
    for (int i = 0; i < radius_index_pairs.size(); ++i) {
        sorted_sphere_table.row(i) = sphere_table.row(radius_index_pairs[i].second);;
    }
*/


    // Convert Table to SpherePack (Physical units)
    Eigen::MatrixX3f centers_phys(sphere_table.rows(), 3);
    Eigen::VectorXf radii_phys(sphere_table.rows());

    for (int i = 0; i < sphere_table.rows(); ++i) {
        Eigen::Vector3f pos_vox = sphere_table.block<1, 3>(i, 0).transpose();
        centers_phys.row(i) = voxel_grid.origin + (pos_vox.array()+0.5f).matrix() * voxel_grid.voxel_size;
        radii_phys(i) = (sphere_table(i, 3)) * voxel_grid.voxel_size;
    }

    // [DEBUG] Verify final sphere positions and radii
    std::cout << "[DEBUG] Final Sphere Sample (Physical Units):" << std::endl;
    for (int i = 0; i < std::min(5, int(sphere_table.rows())); ++i) {
        std::cout <<"x = " << centers_phys(i, 0) << ", y = " << centers_phys(i, 1) << ", z = " << centers_phys(i, 2) << ", r = " << radii_phys(i) << "\n";
    }
    std::cout << std::endl;

    return SpherePack(centers_phys, radii_phys);
}


/**
 * Construct a multisphere representation directly from a triangle mesh.
 * Logically equivalent to the Python version, optimized for C++.
 */

// TODO: Allow for min center and radius to be passed in mesh units (not voxels) and convert internally based on voxel size.
SpherePack multisphere_from_mesh(
    const FastMesh& mesh,
    int div = 100,
    int padding = 2,
    int min_center_distance_vox = 4,
    int max_iter = 1,
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


    // 1. Convert to VoxelGrid (This would call your mesh voxelizer)
    // For parity with your logic, we assume the voxelization results in a grid
    // centered on the mesh centroid with the requested padding.
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
        // 3. Perform Reconstruction
        sp = multisphere_from_voxels(
            voxel_grid,
            min_center_distance_vox,
            max_iter,
            min_radius_vox,
            precision,
            max_spheres,
            show_progress,
            sphere_table_vox // initial spheres in voxel units

        );


    } else {
        std::cout << "No initial sphere table provided. Starting with an empty table." << std::endl;
        sp = multisphere_from_voxels(
            voxel_grid,
            min_center_distance_vox,
            max_iter,
            min_radius_vox,
            precision,
            max_spheres,
            show_progress
        );
    }


    // 4. Boundary Adjustment (Optional)
    if (confine_mesh) {
#ifdef HAVE_MANIFOLD
        // We only include this if the Manifold/Nanoflann logic is enabled
        return utils::adjust_spheres_to_stl_boundary(sp, mesh.vertices);
#else
        constrain_radii_to_sdf(sp, mesh);
        return sp;
#endif
    }

    return sp;
}


#endif