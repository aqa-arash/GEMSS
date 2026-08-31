/**
 * @file main_mesh.cpp
 * @brief Multisphere reconstruction from mesh files using multisphere-cpp.
 *
 * Loads mesh files, reconstructs sphere packs, and exports results.
 *
 * @author Arash Moradian
 * @date 2026-03-09
 */

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

// Project headers
#include "GEMSS/GEMSS-interface.h"

using namespace GEMSS;

void validate_import(const SpherePack& original, const SpherePack& imported, const std::string& format, bool check_physics) {
    if (original.num_spheres() != imported.num_spheres()) {
        std::cerr << "[Fail] " << format << ": Sphere count mismatch.\n";
        return;
    }

    const float tol = 1e-4f;
    bool geo_match = original.centers.isApprox(imported.centers, tol) && 
                     original.radii.isApprox(imported.radii, tol);

    if (!geo_match) {
        std::cerr << "[Fail] " << format << ": Geometric data mismatch.\n";
        return;
    }

    if (check_physics) {
        bool phys_match = 
            std::abs(original.precision - imported.precision) < tol &&
            std::abs(original.mass - imported.mass) < tol &&
            std::abs(original.density - imported.density) < tol &&
            std::abs(original.bounding_radius - imported.bounding_radius) < tol &&
            original.center_of_mass.isApprox(imported.center_of_mass, tol) &&
            original.inertia_tensor.isApprox(imported.inertia_tensor, tol) &&
            original.principal_axes.isApprox(imported.principal_axes, tol) &&
            original.principal_moments.isApprox(imported.principal_moments, tol) &&
            original.is_2d == imported.is_2d;

        if (!phys_match) {
            std::cerr << "[Fail] " << format << ": Physics data mismatch.\n";
            return;
        }
    }

    std::cout << "[Pass] " << format << " import validated successfully.\n";
}

int main() {
    std::cout << "--- Multisphere Reconstruction Test ---" << std::endl;

    std::vector<std::string> models = { "example_mesh.stl" };

    for (const auto& model_name : models) {
        STLMesh example_mesh = load_mesh(model_name);

        GEMSS::MultisphereConfig config;
        config.div = 400; 
        config.padding = 2; 
        config.search_window = 10; 
        config.min_center_distance_rel = 0.5f; 
        config.min_radius_vox = 12; 
        config.precision_target = 0.99f; 
        config.max_spheres = 10000; 
        config.show_progress = true; 
        config.confine_mesh = false; 
        config.initial_sphere_table = Eigen::MatrixXf(0,4); 
        config.compute_physics = 1; 
        config.prune_isolated_spheres = true; 

        SpherePack single_sp = multisphere_from_mesh(
            example_mesh,
            config
        );
        
        export_to_csv(single_sp, model_name + "_recon.csv");
        export_to_vtk(single_sp, model_name + "_recon.vtk");

        print_sphere_pack_info(single_sp);

        auto binary_grid = mesh_to_binary_grid(example_mesh, config); 
        compute_multisphere_physics(single_sp, binary_grid); 
    
        print_sphere_pack_info(single_sp);

        export_sphere_pack_json(single_sp, model_name + "_recon.json");

        SpherePack sp_csv = import_sphere_pack(model_name + "_recon.csv");
        SpherePack sp_vtk = import_sphere_pack(model_name + "_recon.vtk");
        SpherePack sp_json = import_sphere_pack(model_name + "_recon.json");

        validate_import(single_sp, sp_csv, "CSV", false);
        validate_import(single_sp, sp_vtk, "VTK", false);
        validate_import(single_sp, sp_json, "JSON", true);
    }

    return 0;
}