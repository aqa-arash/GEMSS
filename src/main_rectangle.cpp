/**
 * @file main_rectangle.cpp
 * @brief Synthetic rectangular box test for multisphere-cpp reconstruction.
 *
 * Generates a synthetic voxel rectangular box with constant global dimensions,
 * sweeps across multiple voxel sizes while maintaining a constant physical 
 * minimum radius, runs multisphere reconstruction, and exports results.
 *
 * @author Arash Moradian
 * @date 2026-03-24
 */

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include <chrono>
#include "GEMSS/GEMSS-interface.h"


using namespace GEMSS;

int main() {
    std::cout << "--- Multisphere Rectangle Test ---" << std::endl;

    // Define physical dimensions to maintain constant global geometry
    const double length_x = 30.0;
    const double length_y = 20.0;
    const double length_z = 10.0;
    
    // Define physical padding 
    const double pad_x = 0.0;
    const double pad_y = 0.0;
    const double pad_z = 0.0;

    // Define constant physical target for search window and minimum radius
    const float target_min_radius_phys = 0.3f;

    std::vector<double> voxel_sizes = {0.1};

    for (double v_size : voxel_sizes) {
        std::cout << "\n======================================" << std::endl;
        std::cout << "[1/3] Creating grid for v_size = " << v_size << std::endl;

        // Calculate total grid dimensions in voxels
        int nx = static_cast<int>(std::round((length_x + 2.0 * pad_x) / v_size));
        int ny = static_cast<int>(std::round((length_y + 2.0 * pad_y) / v_size));
        int nz = static_cast<int>(std::round((length_z + 2.0 * pad_z) / v_size));

        VoxelGrid<uint8_t> rect(nx, ny, nz, v_size);

        // Calculate start and end indices for the solid rectangle
        int min_x = static_cast<int>(std::round(pad_x / v_size));
        int max_x = static_cast<int>(std::round((pad_x + length_x) / v_size));
        int min_y = static_cast<int>(std::round(pad_y / v_size));
        int max_y = static_cast<int>(std::round((pad_y + length_y) / v_size));
        int min_z = static_cast<int>(std::round(pad_z / v_size));
        int max_z = static_cast<int>(std::round((pad_z + length_z) / v_size));

        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < ny; ++y) {
                for (int z = 0; z < nz; ++z) {
                    if (x >= min_x && x < max_x && y >= min_y && y < max_y && z >= min_z && z < max_z) {
                        rect(x, y, z) = true;
                    }
                }
            }
        }

        std::cout << "[2/3] Running reconstruction..." << std::endl;
        
        // Calculate equivalent voxel counts for the current resolution
        int dynamic_radius_vox = std::max(2, static_cast<int>(std::round(target_min_radius_phys / v_size)));
        
        GEMSS::MultisphereConfig config;
        config.minimum_radius_real = target_min_radius_phys;
        config.min_radius_vox = dynamic_radius_vox;
        config.search_window = dynamic_radius_vox;
        config.precision_target = 0.99f;
        config.max_spheres = 20000; 
        config.show_progress = true;
        config.confine_mesh = false;
        config.compute_physics = 1;

        auto start = std::chrono::high_resolution_clock::now();
        SpherePack rect_sp = multisphere_from_voxels(rect, config);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = end - start; 
        std::cout << " ========== Runtime = " << duration.count()/1000000 << "==========" << std::endl; 

        std::cout << "\nReconstruction Complete for v_size = " << v_size << std::endl;
        print_sphere_pack_info(rect_sp);

        std::cout << "[3/3] Exporting results..." << std::endl;
        
        // Format resolution suffix for output files
        std::string suffix = std::to_string(v_size);
        suffix.erase(suffix.find_last_not_of('0') + 1, std::string::npos);
        if (suffix.back() == '.') {
            suffix += "0";
        }

        std::string csv_file = "reconstructed_rectangle_" + suffix + ".csv";
        std::string stl_file = "original_rectangle_" + suffix + ".stl";
        std::string vtk_file = "reconstructed_rectangle_" + suffix + ".vtk";

        export_to_csv(rect_sp, csv_file);
        save_mesh_to_stl(grid_to_mesh(rect), stl_file);
        export_to_vtk(rect_sp, vtk_file);
    }

    return 0;
}