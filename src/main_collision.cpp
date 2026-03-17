/**
 * @file main_collision.cpp
 * @brief Demonstrates sphere overlap and reconstruction using multisphere-cpp.
 *
 * Generates voxel grids for two overlapping spheres, computes distance transforms,
 * reconstructs sphere packs, and exports results to CSV and VTK files.
 *
 * @author Arash Moradian
 * @date 2026-03-09
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <Eigen/Dense>

// Project headers
#include "multisphere-interface.h"

using namespace MSS;

namespace fs = std::filesystem;

/**
 * @brief Entry point for sphere overlap and reconstruction demo.
 *
 * Generates two spheres with varying center distances, computes distance transforms,
 * reconstructs sphere packs, and exports results.
 *
 * @return int Exit code.
 */
int main() {
    // 1. Directory Setup
    // ------------------
    std::string base_dir = "overlap_study_output";
    std::string csv_dir  = base_dir + "/csv_data";
    std::string bool_dir = base_dir + "/geometries/boolean";
    std::string edt_dir  = base_dir + "/geometries/distance_transform";

    fs::create_directories(csv_dir);
    fs::create_directories(bool_dir);
    fs::create_directories(edt_dir);

    // 2. Simulation Parameters
    // ------------------------
    const int nx = 64, ny = 64, nz = 64;
    const float v_size = 1.0f;
    const float R_large = 20.0f;
    const float R_small = 7.0f;
    const Eigen::Vector3f center_large(R_large+1, 32.0f, 32.0f);

    std::ofstream master_csv(csv_dir + "/summary_results.csv");
    master_csv << "step,center_distance,found_sphere_idx,recon_x,recon_y,recon_z,recon_r\n";

    int step_count = 0;
    for (double dist = R_large+R_small; dist >= R_large-R_small; dist -= 1.0) {
        // 3. Generate Geometry: Union of Two Spheres
        // ------------------------------------------
        VoxelGrid<bool> grid(nx, ny, nz, v_size);
        std::vector<uint8_t> bool_buffer(nx * ny * nz, 0);

        Eigen::Vector3f center_small(center_large.x() + (float)dist, center_large.y(), center_large.z());

        grid.sphere_kernel(center_large.x(), center_large.y(), center_large.z(), R_large, true);
        grid.sphere_kernel(center_small.x(), center_small.y(), center_small.z(), R_small, true);

        // Fill the uint8 buffer for cnpy
        for(size_t i = 0; i < grid.data.size(); ++i) {
            if (grid.data[i]) bool_buffer[i] = 1;
        }

        // 4. Distance Transform & Export
        // -----------------------------
        VoxelGrid<float> edt_grid = grid.distance_transform();


        // 5. Reconstruction

        MSS::MultisphereConfig config;
        config.min_center_distance_vox = 2;
        config.min_radius_vox = 2;
        config.precision_target = 0.99f;
        config.max_spheres = 100000;
        config.show_progress = false;
        // -----------------
        SpherePack pack = multisphere_from_voxels(grid,
                                        config);

        // 6. Export Results
        // -----------------
        for (size_t i = 0; i < pack.num_spheres(); ++i) {
            master_csv << step_count << ","
                       << dist << ","
                       << i << ","
                       << pack.centers(i, 0) << ","
                       << pack.centers(i, 1) << ","
                       << pack.centers(i, 2) << ","
                       << pack.radii(i) << "\n";
        }

        export_to_vtk(pack,
            csv_dir + "/spheres_step_" + std::to_string(step_count) + ".vtk");

        std::cout << "Step " << step_count << " (Dist: " << dist << ") processed." << std::endl;
        step_count++;
    }

    master_csv.close();
    return 0;
}