/**
 * @file main_mesh.cpp
 * @brief Multisphere reconstruction from mesh files using multisphere-cpp.
 *
 * Loads mesh files, reconstructs sphere packs, and exports results.
 * Visualization is enabled if VTK is available.
 *
 * @author Arash Moradian
 * @date 2026-03-09
 */

#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>

// Project headers
#include "multisphere_datatypes.hpp"
#include "multisphere_io.hpp"
#include "multisphere_reconstruction.hpp" // Contains the main algorithm logic

/**
 * @brief Entry point for mesh-based multisphere reconstruction.
 *
 * Loads mesh files, reconstructs sphere packs, and exports/visualizes results.
 *
 * @return int Exit code.
 */
int main() {
    std::cout << "--- Multisphere Reconstruction Test ---" << std::endl;

    // List of mesh files to process
    std::vector<std::string> models = { "xyzrgb_dragon.stl" };

    for (const auto& model_name : models) {
        // 1. Load mesh from file
        FastMesh example_mesh = load_mesh_fast(model_name);

        // --- DEBUG SNIPPET ---
        // Uncomment to save mesh for debugging
        // save_mesh_to_stl(example_mesh, "debug_output.stl");
        // std::cout << "[DEBUG] Check 'debug_output.stl' against 'example_mesh.stl' now." << std::endl;

        // 2. Run the Reconstruction Algorithm
        std::cout << "[2/3] Running reconstruction..." << std::endl;

        SpherePack single_sp = multisphere_from_mesh(
            example_mesh,
            400,    // div
            2,      // padding
            10,     // min_center_distance_vox
            10,     // max_iter
            4,      // min_radius_vox
            0.99,   // precision_target
            100000, // max_spheres
            true,   // show_progress
            false   // confine_mesh
        );

        std::cout << "\nReconstruction Complete!" << std::endl;
        std::cout << "--Single Sphere : \n Spheres found: " << single_sp.num_spheres() << std::endl;
        std::cout << "Max radius: " << single_sp.max_radius() << " units" << std::endl;

        // 3. Visualization or Export

        std::cout << "[3/3] VTK not enabled. Skipping visualization." << std::endl;
        export_to_csv(single_sp, model_name + ".csv");
        export_to_vtk(single_sp, model_name + ".vtk");
        std::cout << "To see the result, export to CSV or recompile with VTK enabled." << std::endl;
    }

    return 0;
}