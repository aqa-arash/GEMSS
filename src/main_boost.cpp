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
    std::string base_dir = "Boost_study_output";

    //////////////////////////////////////////////////////////////////////////////////////
    // case one : Ellipsoid with two large spheres at the center with high overlap to kickstart the reconstruction
    /////////////////////////////////////////////////////////////////////////////////////

    std::string case_dir = base_dir + "/case_ellipsoid_overlap";
    std::string csv_dir  = case_dir + "/csv_data";
    std::string bool_dir = case_dir + "/geometries/boolean";
    std::string edt_dir  = case_dir + "/geometries/distance_transform";

    fs::create_directories(csv_dir);
    fs::create_directories(bool_dir);
    fs::create_directories(edt_dir);

    // 2. Simulation Parameters
    const int nx = 256, ny = 256, nz = 256;
    Eigen::Vector3f center_point = Eigen::Vector3f(nx/2.0f, ny/2.0f, nz/2.0f);
    const float v_size = 1.0f; // voxel size
       // Initialize Boolean Grid
    VoxelGrid<bool> grid(nx, ny, nz, v_size);
    std::vector<uint8_t> bool_buffer(nx * ny * nz, 0);

    VoxelGrid<float> edt_grid = grid.distance_transform();

    
    const float a = 120.0f;
    const float b = 60.0f;
    const float c = 60.0f;
    // make an ellipsoid with radii a,b,c and then scale it to fit in the grid
 
    
    std::ofstream master_csv(csv_dir + "/summary_results.csv");
    master_csv << "step,center_distance,found_sphere_idx,recon_x,recon_y,recon_z,recon_r\n";




    // 3. Generate Ellipsoid Geometry DO NOT PARALLELIZE
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                float dx = (x - center_point.x()) * v_size;
                float dy = (y - center_point.y()) * v_size;
                float dz = (z - center_point.z()) * v_size;
                float val = (dx*dx)/(a*a) + (dy*dy)/(b*b) + (dz*dz)/(c*c);
                if (val <= 1.0f) {
                    grid(x, y, z) = true;
                }
                else {
                    grid(x, y, z) = false;
                }
            }
        }
    }


    // Fill the uint8 buffer for cnpy
    for(size_t i = 0; i < grid.data.size(); ++i) {
        if (grid.data[i]) bool_buffer[i] = 1;
    }

    // 3. Distance Transform & Export
    VoxelGrid<float> edt_grid = grid.distance_transform(); 

    cnpy::npy_save(bool_dir + "/bool_step_"  + ".npy", 
                    bool_buffer.data(), {(size_t)nx, (size_t)ny, (size_t)nz}, "w");

    cnpy::npy_save(edt_dir + "/edt_step_" + ".npy", 
                    edt_grid.data.data(), {(size_t)nx, (size_t)ny, (size_t)nz}, "w");

    Eigen::MatrixX4f sphere_table = Eigen::MatrixX4f(2, 4); // Empty table for initial reconstruction
    // Add two large spheres close to the center to kickstart the reconstruction

    float center_distance = 20.0f; // Distance from center for initial spheres

    float radius0 = edt_grid((int)(center_point.x() - center_distance), (int)center_point.y(), (int)center_point.z());
    float radius1 = edt_grid((int)(center_point.x() + center_distance), (int)center_point.y(), (int)center_point.z());

    std::cout << "Initial Sphere 0: Center = (" << center_point.x() - center_distance << ", " << center_point.y() << ", " << center_point.z() << "), Radius = " << radius0 << std::endl;
    std::cout << "Initial Sphere 1: Center = (" << center_point.x() + center_distance << ", " << center_point.y() << ", " << center_point.z() << "), Radius = " << radius1 << std::endl;

    sphere_table.row(0) = Eigen::Vector4f(center_point.x() - center_distance, center_point.y(), center_point.z(), radius0);
    sphere_table.row(1) = Eigen::Vector4f(center_point.x() + center_distance, center_point.y(), center_point.z(), radius1);


    // 4. Reconstruction
    SpherePack pack = multisphere_from_voxels(grid,
                                    10, // min_center_distance_vox
                                    100, // max iter
                                    20, // min_radius_vox
                                    0.99, // precision_target
                                    100000, // max_spheres
                                    false,
                                    sphere_table // initial sphere table
                                    );



    // 5. FIXED: Accessing Sphere Data via Eigen Matrices
    for (size_t i = 0; i < pack.num_spheres(); ++i) {
        // pack.centers is MatrixX3f (Rows = spheres, Cols = X,Y,Z)
        // pack.radii is VectorXf
        master_csv  << "final," 
                    << center_distance << "," 
                    << i << "," 
                    << pack.centers(i, 0) << "," // X
                    << pack.centers(i, 1) << "," // Y
                    << pack.centers(i, 2) << "," // Z
                    << pack.radii(i) << "\n";    // Radius
    }

    export_to_vtk(pack, 
        csv_dir + "/spheres_step_" +  ".vtk");

    std::cout  << " (Dist: " << center_distance << ") processed." << std::endl;
    master_csv.close();

    /////////////////////////////////////////////////////////////////////////////////
    // case two : Ellipsoid with large min center distance and initial spheres with minimum overlap
    /////////////////////////////////////////////////////////////////////////////////


    case_dir = base_dir + "/case_ellipsoid_min_dist";
    csv_dir  = case_dir + "/csv_data";
    bool_dir = case_dir + "/geometries/boolean";
    edt_dir  = case_dir + "/geometries/distance_transform";

    fs::create_directories(csv_dir);
    fs::create_directories(bool_dir);
    fs::create_directories(edt_dir);


    std::ofstream master_csv1(csv_dir + "/summary_results.csv");
    master_csv1 << "step,center_distance,found_sphere_idx,recon_x,recon_y,recon_z,recon_r\n";

    // 3. Generate Ellipsoid Geometry (DO NOT PARALLELIZE)
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                float dx = (x - center_point.x()) * v_size;
                float dy = (y - center_point.y()) * v_size;
                float dz = (z - center_point.z()) * v_size;
                float val = (dx*dx)/(a*a) + (dy*dy)/(b*b) + (dz*dz)/(c*c);
                if (val <= 1.0f) {
                    grid(x, y, z) = true;
                }
                else {
                    grid(x, y, z) = false;
                }
            }
        }
    }


    // Fill the uint8 buffer for cnpy
    for(size_t i = 0; i < grid.data.size(); ++i) {
        if (grid.data[i]) bool_buffer[i] = 1;
    }

    // 3. Distance Transform & Export
    edt_grid = grid.distance_transform(); 

    cnpy::npy_save(bool_dir + "/bool_step_"  + ".npy", 
                    bool_buffer.data(), {(size_t)nx, (size_t)ny, (size_t)nz}, "w");

    cnpy::npy_save(edt_dir + "/edt_step_" + ".npy", 
                    edt_grid.data.data(), {(size_t)nx, (size_t)ny, (size_t)nz}, "w");

    sphere_table = Eigen::MatrixX4f(2, 4); // Empty table for initial reconstruction
    // Add two large spheres close to the center to kickstart the reconstruction

    center_distance = 20.0f; // Distance from center for initial spheres
    
    radius0 = edt_grid((int)(center_point.x() - center_distance), (int)center_point.y(), (int)center_point.z());
    radius1 = edt_grid((int)(center_point.x() + center_distance), (int)center_point.y(), (int)center_point.z());

    float temp_center_distance = center_distance;
    float temp_r0 = radius0;
    float temp_r1 = radius1;
    bool overlap_exist = radius0 + radius1 - 2 * center_distance > 0.0f;
    for (int i = 0; overlap_exist  == true; ++i) {
        center_distance = temp_center_distance;
        radius0 = temp_r0;
        radius1 = temp_r1;
        temp_center_distance += 3.0f; // Increase distance until nearly no overlap
        temp_r0 = edt_grid((int)(center_point.x() - temp_center_distance), (int)center_point.y(), (int)center_point.z());
        temp_r1 = edt_grid((int)(center_point.x() + temp_center_distance), (int)center_point.y(), (int)center_point.z());
        overlap_exist = temp_r0 + temp_r1 - 2 * temp_center_distance > 0.0f;
    }



    std::cout << "Initial Sphere 0: Center = (" << center_point.x() - center_distance << ", " << center_point.y() << ", " << center_point.z() << "), Radius = " << radius0 << std::endl;
    std::cout << "Initial Sphere 1: Center = (" << center_point.x() + center_distance << ", " << center_point.y() << ", " << center_point.z() << "), Radius = " << radius1 << std::endl;

    sphere_table.row(0) = Eigen::Vector4f(center_point.x() - center_distance, center_point.y(), center_point.z(), radius0);
    sphere_table.row(1) = Eigen::Vector4f(center_point.x() + center_distance, center_point.y(), center_point.z(), radius1);


    // 4. Reconstruction
    SpherePack pack2 = multisphere_from_voxels(grid,
                                    center_distance + 2, // min_center_distance_vox
                                    20, // max iter
                                    20, // min_radius_vox
                                    0.99, // precision_target
                                    100000, // max_spheres
                                    false,
                                    sphere_table // initial sphere table
                                    );



    // 5. FIXED: Accessing Sphere Data via Eigen Matrices
    for (size_t i = 0; i < pack2.num_spheres(); ++i) {
        // pack2.centers is MatrixX3f (Rows = spheres, Cols = X,Y,Z)
        // pack2.radii is VectorXf
        master_csv1  << "final," 
                    << center_distance << "," 
                    << i << "," 
                    << pack2.centers(i, 0) << "," // X
                    << pack2.centers(i, 1) << "," // Y
                    << pack2.centers(i, 2) << "," // Z
                    << pack2.radii(i) << "\n";    // Radius
    }

    export_to_vtk(pack2, 
        csv_dir + "/spheres_step_" +  ".vtk");

    std::cout  << " (Dist: " << center_distance << ") processed." << std::endl;
    master_csv1.close();


    ///////////////////////////////////////////////////////////////////////////////
    // case three : cube 
    ///////////////////////////////////////////////////////////////////////////////

    
    case_dir = base_dir + "/case_cube";
    csv_dir  = case_dir + "/csv_data";
    bool_dir = case_dir + "/geometries/boolean";
    edt_dir  = case_dir + "/geometries/distance_transform";

    fs::create_directories(csv_dir);
    fs::create_directories(bool_dir);
    fs::create_directories(edt_dir);

    const float cube_size = 200.0f;
    const float half_s = cube_size / 2.0f;

    
    std::ofstream master_csv2(csv_dir + "/summary_results.csv");
    master_csv2 << "step,center_distance,found_sphere_idx,recon_x,recon_y,recon_z,recon_r\n";
    VoxelGrid<float> medial_axis_grid(nx, ny, nz, v_size);
    
 // DO NOT PARALLELIZE (FAILS ON BOOLEAN)
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                float dx = std::abs(x - center_point.x()) * v_size;
                float dy = std::abs(y - center_point.y()) * v_size;
                float dz = std::abs(z - center_point.z()) * v_size;


                if (dx <= half_s && dy <= half_s && dz <= half_s) {
                    grid(x, y, z) = true;
                    
                    // Medial Axis Logic
                    // 1. Calculate distance to nearest faces
                    float dist_x = half_s - std::abs(x - center_point.x());
                    float dist_y = half_s - std::abs(y - center_point.y());
                    float dist_z = half_s - std::abs(z - center_point.z());
                    // 2. The medial axis exists where the two SMALLEST distances are nearly equal
                    // We sort them or check pairs:
                    if (dist_x >1 && dist_y > 1 && dist_z > 1) { // Only consider points inside the cube for medial axis
                        bool on_axis = false;
                        float eps = v_size * 0.8f; // Adaptive threshold based on resolution
                        if (std::abs(dist_x - dist_y) < eps && dist_x <= dist_z) on_axis = true;
                        else if (std::abs(dist_x - dist_z) < eps && dist_x <= dist_y) on_axis = true;
                        else if (std::abs(dist_y - dist_z) < eps && dist_y <= dist_x) on_axis = true;
                        if (on_axis) {
                            medial_axis_grid(x, y, z) = 1.0f; // Store distance to surface on medial axis
                        }
                        else {
                            medial_axis_grid(x, y, z) = 0.0f; // Not part of the medial axis
                        }
                    }
                    
                }
                else {
                    grid(x, y, z) = false;
                }
            }
        }
    }

    export_voxel_grid_to_vtk(medial_axis_grid, csv_dir + "/medial_axis.vtk");
    FastMesh Medial = grid_to_mesh(medial_axis_grid);
    save_mesh_to_stl(Medial, csv_dir + "/surface_mesh.stl");


    // Fill the uint8 buffer for cnpy
    #pragma omp parallel for 
    for(size_t i = 0; i < grid.data.size(); ++i) {
        if (grid.data[i]) bool_buffer[i] = 1;
    }   

    // 3. Distance Transform & Export
    edt_grid = grid.distance_transform(); 

    

    cnpy::npy_save(bool_dir + "/bool_step_"  + ".npy", 
                    bool_buffer.data(), {(size_t)nx, (size_t)ny, (size_t)nz}, "w");

    cnpy::npy_save(edt_dir + "/edt_step_" + ".npy", 
                    edt_grid.data.data(), {(size_t)nx, (size_t)ny, (size_t)nz}, "w");

    
    
    // 4. Reconstruction
    SpherePack pack3 = multisphere_from_voxels(grid,
                                    10, // min_center_distance_vox
                                    10, // max iter
                                    20, // min_radius_vox
                                    0.99, // precision_target
                                    100000, // max_spheres
                                    false
                                    );



    // 5. FIXED: Accessing Sphere Data via Eigen Matrices
    for (size_t i = 0; i < pack3.num_spheres(); ++i) {
        // pack3.centers is MatrixX3f (Rows = spheres, Cols = X,Y,Z)
        // pack3.radii is VectorXf
        master_csv2  << "final," 
                    << i << "," 
                    << pack3.centers(i, 0) << "," // X
                    << pack3.centers(i, 1) << "," // Y
                    << pack3.centers(i, 2) << "," // Z
                    << pack3.radii(i) << "\n";    // Radius
    }

    export_to_vtk(pack3, 
        csv_dir + "/spheres_step_" +  ".vtk");

    std::cout  << ") processed." << std::endl;
    master_csv2.close();

    
    return 0;
}