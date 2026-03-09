#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem> 
#include <Eigen/Dense>
#include "multisphere_datatypes.hpp"
#include "multisphere_reconstruction.hpp"
#include "multisphere_visualization.hpp"
#include "cnpy.h" 

namespace fs = std::filesystem;

int main() {
    // 1. Directory Setup
    std::string base_dir = "overlap_study_output";
    std::string csv_dir  = base_dir + "/csv_data";
    std::string bool_dir = base_dir + "/geometries/boolean";
    std::string edt_dir  = base_dir + "/geometries/distance_transform";

    fs::create_directories(csv_dir);
    fs::create_directories(bool_dir);
    fs::create_directories(edt_dir);

    // 2. Simulation Parameters
    const int nx = 64, ny = 64, nz = 64;
    const float v_size = 1.0f;
    const float R_large = 20.0f;
    const float R_small = 7.0f;
    const Eigen::Vector3f center_large(R_large+1, 32.0f, 32.0f);

    std::ofstream master_csv(csv_dir + "/summary_results.csv");
    master_csv << "step,center_distance,found_sphere_idx,recon_x,recon_y,recon_z,recon_r\n";

    int step_count = 0;
    for (double dist = R_large+R_small; dist >= R_large-R_small; dist -= 1.0) {
        // Initialize Boolean Grid
        VoxelGrid<bool> grid(nx, ny, nz, v_size);
        std::vector<uint8_t> bool_buffer(nx * ny * nz, 0);

        // Generate geometry: Union of two spheres
        Eigen::Vector3f center_small(center_large.x() + (float)dist, center_large.y(), center_large.z());
        
        // Use the built-in sphere_kernel for efficiency
        grid.sphere_kernel(center_large.x(), center_large.y(), center_large.z(), R_large, true);
        grid.sphere_kernel(center_small.x(), center_small.y(), center_small.z(), R_small, true);

        // Fill the uint8 buffer for cnpy
        for(size_t i = 0; i < grid.data.size(); ++i) {
            if (grid.data[i]) bool_buffer[i] = 1;
        }

        // 3. Distance Transform & Export
        VoxelGrid<float> edt_grid = grid.distance_transform(); 

        cnpy::npy_save(bool_dir + "/bool_step_" + std::to_string(step_count) + ".npy", 
                       bool_buffer.data(), {(size_t)nx, (size_t)ny, (size_t)nz}, "w");

        cnpy::npy_save(edt_dir + "/edt_step_" + std::to_string(step_count) + ".npy", 
                       edt_grid.data.data(), {(size_t)nx, (size_t)ny, (size_t)nz}, "w");

        // 4. Reconstruction
        SpherePack pack = multisphere_from_voxels(grid,
                                        2, // min_center_distance_vox
                                        5, // max iter
                                        2, // min_radius_vox
                                        0.99, // precision_target
                                        100000, // max_spheres
                                           false);



        // 5. FIXED: Accessing Sphere Data via Eigen Matrices
        for (size_t i = 0; i < pack.num_spheres(); ++i) {
            // pack.centers is MatrixX3f (Rows = spheres, Cols = X,Y,Z)
            // pack.radii is VectorXf
            master_csv << step_count << "," 
                       << dist << "," 
                       << i << "," 
                       << pack.centers(i, 0) << "," // X
                       << pack.centers(i, 1) << "," // Y
                       << pack.centers(i, 2) << "," // Z
                       << pack.radii(i) << "\n";    // Radius
        }

        export_to_vtk(pack, 
            csv_dir + "/spheres_step_" + std::to_string(step_count) + ".vtk");

        std::cout << "Step " << step_count << " (Dist: " << dist << ") processed." << std::endl;
        step_count++;
    }

    master_csv.close();
    return 0;
}