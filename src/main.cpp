/**
 * @file main.cpp
 * @brief Synthetic geometry tests for multisphere-cpp reconstruction.
 *
 * Generates various synthetic voxel geometries (sphere, double sphere, rectangle, L-shape, torus, cluster, sputnik, tumor, cheese, ice cream cone),
 * runs multisphere reconstruction, and exports results.
 * Visualization is enabled if VTK is available.
 *
 * @author Arash Moradian
 * @date 2026-03-09
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "multisphere-interface.h"

using namespace MSS;


void print_sphere_pack_info(const SpherePack& sp) {
    std::cout << "Sphere Pack Info:" << std::endl;
    std::cout << "      Number of spheres: " << sp.num_spheres() << std::endl;
    std::cout << "      Max radius: " << sp.max_radius() << " units" << std::endl;
    std::cout << "      Min radius: " << sp.min_radius() << " units" << std::endl;
    std::cout << "      Volume of union: " << sp.volume << " units^3" << std::endl;
    std::cout << "      Center of mass: " << sp.center_of_mass.transpose() << " units" << std::endl;
    std::cout << "      Principal moments: " << sp.principal_moments.transpose() << " units^5" << std::endl;
    std::cout << "      Principal axes:\n" << sp.principal_axes << std::endl;
}


/**
 * @brief Entry point for synthetic geometry reconstruction demo.
 *
 * Generates multiple synthetic voxel grids, reconstructs sphere packs, and exports/visualizes results.
 *
 * @return int Exit code.
 */
int main() {
    std::cout << "--- Multisphere Synthetic Geometry Test ---" << std::endl;

    // 1. Create synthetic voxel grids for various shapes
    // --------------------------------------------------
    double v_size = 1.0;
    int nx = 64, ny = 64, nz = 64;

    // Single Sphere
    VoxelGrid<bool> single_sphere(nx, ny, nz, v_size);
    int cx = 32, cy = 32, cz = 32, r = 15;
    for (int x = 0; x < nx; ++x)
        for (int y = 0; y < ny; ++y)
            for (int z = 0; z < nz; ++z)
                if (std::sqrt(std::pow(x-cx,2)+std::pow(y-cy,2)+std::pow(z-cz,2)) <= r)
                    single_sphere.data[x * ny * nz + y * nz + z] = true;

    // Double Sphere
    VoxelGrid<bool> double_sphere(nx, ny, nz, v_size);
    struct Ball { double x, y, z, r; };
    std::vector<Ball> test_geometry = {
        {30.0, 32.0, 32.0, 14.0},
        {35.0, 32.0, 32.0, 14.0},
        {40.0, 32.0, 32.0, 14.0}
    };
    for (int x = 0; x < nx; ++x)
        for (int y = 0; y < ny; ++y)
            for (int z = 0; z < nz; ++z) {
                bool is_inside = false;
                for (const auto& b : test_geometry)
                    if (std::pow(x-b.x,2)+std::pow(y-b.y,2)+std::pow(z-b.z,2) <= b.r*b.r) {
                        is_inside = true;
                        break;
                    }
                if (is_inside)
                    double_sphere.data[x * ny * nz + y * nz + z] = true;
            }

    // Rectangle
    VoxelGrid<bool> rectangle(nx, ny, nz, v_size);
    int min_x = 12, max_x = 55, min_y = 22, max_y = 42, min_z = 22, max_z = 42;
    for (int x = 0; x < nx; ++x)
        for (int y = 0; y < ny; ++y)
            for (int z = 0; z < nz; ++z)
                if (x >= min_x && x < max_x && y >= min_y && y < max_y && z >= min_z && z < max_z)
                    rectangle.data[x * ny * nz + y * nz + z] = true;

    // L-Shape
    VoxelGrid<bool> l_shape(nx, ny, nz, v_size);
    for (int x = 0; x < nx; ++x)
        for (int y = 0; y < ny; ++y)
            for (int z = 0; z < nz; ++z) {
                bool vert = (x >= 20 && x < 30) && (y >= 20 && y < 50) && (z >= 20 && z < 30);
                bool horz = (x >= 20 && x < 50) && (y >= 20 && y < 30) && (z >= 20 && z < 30);
                if (vert || horz)
                    l_shape.data[x * ny * nz + y * nz + z] = true;
            }

    // Torus
    VoxelGrid<bool> torus(nx, ny, nz, v_size);
    double major_R = 20.0, minor_r = 6.0;
    for (int x = 0; x < nx; ++x)
        for (int y = 0; y < ny; ++y)
            for (int z = 0; z < nz; ++z) {
                double dx = x - 32.0, dy = y - 32.0, dz = z - 32.0;
                double len_xy = std::sqrt(dx*dx + dy*dy);
                double tube_dist = std::sqrt(std::pow(len_xy - major_R, 2) + dz*dz);
                if (tube_dist <= minor_r)
                    torus.data[x * ny * nz + y * nz + z] = true;
            }

    // Phantom Peak Cluster
    VoxelGrid<bool> cluster(nx, ny, nz, v_size);
    struct BallData { double x, y, z, r; };
    std::vector<BallData> cluster_spheres = {
        {25.0, 25.0, 32.0, 18.0}, {25.0, 39.0, 32.0, 18.0},
        {39.0, 25.0, 32.0, 18.0}, {39.0, 39.0, 32.0, 18.0}
    };
    for (int x = 0; x < nx; ++x)
        for (int y = 0; y < ny; ++y)
            for (int z = 0; z < nz; ++z) {
                bool is_inside = false;
                for (const auto& b : cluster_spheres)
                    if (std::pow(x-b.x,2)+std::pow(y-b.y,2)+std::pow(z-b.z,2) <= b.r*b.r) {
                        is_inside = true;
                        break;
                    }
                if (is_inside)
                    cluster(x, y, z) = true;
            }

    // Sputnik
    VoxelGrid<bool> sputnik(nx, ny, nz, v_size);
    struct SputnikBall { double x, y, z, r; };
    std::vector<SputnikBall> sputnik_geo = {
        {32.0, 32.0, 32.0, 20.0},
        {32.0+24.0, 32.0, 32.0, 8.0}, {32.0-24.0, 32.0, 32.0, 8.0},
        {32.0, 32.0+24.0, 32.0, 8.0}, {32.0, 32.0-24.0, 32.0, 8.0},
        {32.0, 32.0, 32.0+24.0, 8.0}, {32.0, 32.0, 32.0-24.0, 8.0}
    };
    for (int x = 0; x < nx; ++x)
        for (int y = 0; y < ny; ++y)
            for (int z = 0; z < nz; ++z) {
                bool is_inside = false;
                for (const auto& b : sputnik_geo)
                    if (std::pow(x-b.x,2)+std::pow(y-b.y,2)+std::pow(z-b.z,2) <= b.r*b.r) {
                        is_inside = true;
                        break;
                    }
                if (is_inside)
                    sputnik(x, y, z) = true;
            }

    // Tumor
    VoxelGrid<bool> tumor(nx, ny, nz, v_size);
    struct TumorBall { double x, y, z, r; };
    std::vector<TumorBall> tumor_geo = {
        {32.0, 32.0, 32.0, 20.0}, {44.0, 32.0, 32.0, 15.0}
    };
    for (int x = 0; x < nx; ++x)
        for (int y = 0; y < ny; ++y)
            for (int z = 0; z < nz; ++z) {
                bool is_inside = false;
                for (const auto& b : tumor_geo)
                    if (std::pow(x-b.x,2)+std::pow(y-b.y,2)+std::pow(z-b.z,2) <= b.r*b.r) {
                        is_inside = true;
                        break;
                    }
                if (is_inside)
                    tumor(x, y, z) = true;
            }

    // Swiss Cheese
    VoxelGrid<bool> cheese(nx, ny, nz, v_size);
    int box_min = 12, box_max = 52;
    struct Hole { double x, y, z, r; };
    std::vector<Hole> holes = {
        {20.0, 20.0, 52.0, 12.0}, {45.0, 45.0, 12.0, 15.0}, {32.0, 32.0, 32.0, 10.0}
    };
    for (int x = 0; x < nx; ++x)
        for (int y = 0; y < ny; ++y)
            for (int z = 0; z < nz; ++z) {
                bool inside_box = (x >= box_min && x < box_max) &&
                                  (y >= box_min && y < box_max) &&
                                  (z >= box_min && z < box_max);
                if (!inside_box) continue;
                bool inside_hole = false;
                for (const auto& h : holes)
                    if (std::pow(x-h.x,2)+std::pow(y-h.y,2)+std::pow(z-h.z,2) <= h.r*h.r) {
                        inside_hole = true;
                        break;
                    }
                if (!inside_hole)
                    cheese(x, y, z) = true;
            }
    save_mesh_to_stl(grid_to_mesh(cheese), "swiss_cheese.stl");

    // Ice Cream Cone
    VoxelGrid<uint8_t> ice_cream(nx, ny, nz, v_size);
    double cone_tip_z = 2.0, cone_base_z = 50.0, cone_base_r = 10.0, scoop_r = 12.0;
    for (int x = 0; x < nx; ++x)
        for (int y = 0; y < ny; ++y)
            for (int z = 0; z < nz; ++z) {
                double dx = x - 32.0, dy = y - 32.0;
                bool in_cone = false;
                if (z >= cone_tip_z && z <= cone_base_z) {
                    double height_fraction = (z - cone_tip_z) / (cone_base_z - cone_tip_z);
                    double current_r = height_fraction * cone_base_r;
                    if ((dx*dx + dy*dy) <= (current_r * current_r))
                        in_cone = true;
                }
                double dz_sphere = z - cone_base_z;
                bool in_sphere = (dx*dx + dy*dy + dz_sphere*dz_sphere) <= (scoop_r * scoop_r);
                if (in_cone || in_sphere)
                    ice_cream(x, y, z) = true;
            }

            save_mesh_to_stl(grid_to_mesh(ice_cream), "ice_cream_cone.stl");

    // 2. Run multisphere reconstruction for each geometry
    // ---------------------------------------------------
    std::cout << "[2/3] Running reconstruction..." << std::endl;

    MSS::MultisphereConfig config;
        config.min_center_distance_vox = 4;
        config.min_radius_vox = 2;
        config.precision_target = 0.95f;
        config.max_spheres = 10;
        config.show_progress = true;
        config.confine_mesh = true;
        config.compute_physics = true;

    // voxelgrid, min_center_distance_vox, min_radius_vox, precision_target, max_spheres, use_boost, initial_sphere_table
    SpherePack single_sp = multisphere_from_voxels(single_sphere, config);
    SpherePack double_sp = multisphere_from_voxels(double_sphere, config);
    SpherePack rectangle_sp = multisphere_from_voxels(rectangle, config);
    SpherePack l_shape_sp = multisphere_from_voxels(l_shape, config);
    SpherePack torus_sp = multisphere_from_voxels(torus, config);
    SpherePack cluster_sp = multisphere_from_voxels(cluster, config);
    SpherePack sputnik_sp = multisphere_from_voxels(sputnik, config);
    SpherePack tumor_sp = multisphere_from_voxels(tumor, config);
    SpherePack cheese_sp = multisphere_from_voxels(cheese, config);
    SpherePack ice_cream_boost = multisphere_from_voxels(ice_cream, config);

    std::cout << "\nReconstruction Complete!" << std::endl;
    std::cout << "--Single Sphere : " << std::endl;
    print_sphere_pack_info(single_sp);
    std::cout << "--Double Sphere : " << std::endl;
    print_sphere_pack_info(double_sp);
    std::cout << "--Rectangle : " << std::endl;
    print_sphere_pack_info(rectangle_sp);
    std::cout << "--L-Shape : " << std::endl;
    print_sphere_pack_info(l_shape_sp);
    std::cout << "--Torus : " << std::endl;
    print_sphere_pack_info(torus_sp);
    std::cout << "--Cluster : " << std::endl;
    print_sphere_pack_info(cluster_sp);
    std::cout << "--Sputnik : " << std::endl;
    print_sphere_pack_info(sputnik_sp);
    std::cout << "--Tumor : " << std::endl;
    print_sphere_pack_info(tumor_sp);
    std::cout << "--Cheese : " << std::endl;
    print_sphere_pack_info(cheese_sp);
    std::cout << "--Ice Cream Cone (Boost) : " << std::endl;
    print_sphere_pack_info(ice_cream_boost);   


    // 3. Visualization or Export
    // --------------------------
    export_to_csv(single_sp, "reconstructed_spheres.csv");
    export_to_vtk(single_sp, "reconstructed_spheres.vtk");
    export_to_csv(double_sp, "reconstructed_double_spheres.csv");
    export_to_vtk(double_sp, "reconstructed_double_spheres.vtk");
    export_to_csv(rectangle_sp, "reconstructed_rectangle_spheres.csv");
    export_to_vtk(rectangle_sp, "reconstructed_rectangle_spheres.vtk");
    export_to_csv(l_shape_sp, "reconstructed_l_shape.csv");
    export_to_vtk(l_shape_sp, "reconstructed_l_shape.vtk");
    export_to_csv(torus_sp, "reconstructed_torus.csv");
    export_to_vtk(torus_sp, "reconstructed_torus.vtk");
    export_to_csv(cluster_sp, "reconstructed_cluster.csv");
    export_to_vtk(cluster_sp, "reconstructed_cluster.vtk");
    export_to_csv(sputnik_sp, "reconstructed_sputnik.csv");
    export_to_vtk(sputnik_sp, "reconstructed_sputnik.vtk");
    export_to_csv(tumor_sp, "reconstructed_tumor.csv");
    export_to_vtk(tumor_sp, "reconstructed_tumor.vtk");
    export_to_csv(cheese_sp, "reconstructed_cheese.csv");
    export_to_vtk(cheese_sp, "reconstructed_cheese.vtk");
    export_to_csv(ice_cream_boost, "reconstructed_ice_cream_boost.csv");
    export_to_vtk(ice_cream_boost, "reconstructed_ice_cream_boost.vtk");
    std::cout << "To see the result, export to CSV or recompile with VTK enabled." << std::endl;

    return 0;
}