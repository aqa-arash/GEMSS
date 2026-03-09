#include <iostream>
#include <Eigen/Dense>
#include "multisphere_datatypes.hpp"
#include "multisphere_reconstruction.hpp" // Contains the main algorithm logic
#include "multisphere_visualization.hpp"

int main() {
    std::cout << "--- Multisphere Reconstruction Test ---" << std::endl;

    // 1. Create a synthetic VoxelGrid (64x64x64)
    // We'll place a manual "blob" in the center to reconstruct
    double v_size = 1.0;
    int nx = 64, ny = 64, nz = 64;
    VoxelGrid<bool> single_sphere(nx, ny, nz, v_size);
    
    int cx = 32, cy = 32, cz = 32, r = 15;
    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int z = 0; z < nz; ++z) {
                double dist = std::sqrt(std::pow(x-cx, 2) + std::pow(y-cy, 2) + std::pow(z-cz, 2));
                if (dist <= r) {
                    single_sphere.data[x * (64 * 64) + y * 64 + z] = true;
                }
            }
        }
    }

    // Create Geometry: Union of Two Spheres
    // Sphere A: Center (20, 32, 32), Radius 14
    // Sphere B: Center (44, 32, 32), Radius 14
    // Distance 24, Overlap 4 units.
    VoxelGrid<bool> double_sphere(nx, ny, nz, v_size);


    struct Ball { double x, y, z, r; };
    std::vector<Ball> test_geometry = {
        {30.0, 32.0, 32.0, 14.0},
        {35.0, 32.0, 32.0, 14.0},
        {40.0, 32.0, 32.0, 14.0}
    };

    std::cout << "[1/3] Generating geometry..." << std::endl;

    // Fill grid (Union operation)
    // Using serial loop to safely write to vector<bool> proxy
    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int z = 0; z < nz; ++z) {
                bool is_inside = false;
                for (const auto& b : test_geometry) {
                    double dist_sq = std::pow(x - b.x, 2) + 
                                     std::pow(y - b.y, 2) + 
                                     std::pow(z - b.z, 2);
                    if (dist_sq <= b.r * b.r) {
                        is_inside = true;
                        break;
                    }
                }
                
                if (is_inside) {
                    // Use operator() if available, or manual index
                    double_sphere.data[x * ny * nz + y * nz + z] = true; 
                }
            }
        }
    }

    // 2. Create Geometry: Rectangular Cuboid
    // Dimensions: 40 (x) by 20 (y) by 20 (z)
    // The "Plateau" is a line along the X-axis.
    VoxelGrid<bool> rectangle(nx, ny, nz, v_size);    

    
    int min_x = 12, max_x = 55; // Length 40
    int min_y = 22, max_y = 42; // Width 20
    int min_z = 22, max_z = 42; // Height 20

    // Fill the grid
    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int z = 0; z < nz; ++z) {
                // Check bounds
                if (x >= min_x && x < max_x &&
                    y >= min_y && y < max_y &&
                    z >= min_z && z < max_z) {
                    
                    // Set voxel to true (Object)
                    rectangle.data[x * ny * nz + y * nz + z] = true; 
                }
            }
        }
    }

    // 3. Create Geometry: L-Shape (Concave Logic Test)
    // Vertical Bar:   [20..30] x [20..50] x [20..30]
    // Horizontal Bar: [20..50] x [20..30] x [20..30]
    // Common Corner:  [20..30] x [20..30] x [20..30]
    VoxelGrid<bool> l_shape(nx, ny, nz, v_size);
    
    std::cout << "[1/3] Generating L-Shape..." << std::endl;
    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int z = 0; z < nz; ++z) {
                // Vertical Leg
                bool vert = (x >= 20 && x < 30) && (y >= 20 && y < 50) && (z >= 20 && z < 30);
                // Horizontal Leg
                bool horz = (x >= 20 && x < 50) && (y >= 20 && y < 30) && (z >= 20 && z < 30);
                
                if (vert || horz) {
                    l_shape.data[x * ny * nz + y * nz + z] = true;
                }
            }
        }
    }

    // 4. Create Geometry: Torus (Donut)
    // Center (32,32,32). Major Radius R=20, Minor Radius r=6.
    VoxelGrid<bool> torus(nx, ny, nz, v_size);
    double major_R = 20.0;
    double minor_r = 6.0;

    std::cout << "[1/3] Generating Torus..." << std::endl;
    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int z = 0; z < nz; ++z) {
                double dx = x - 32.0;
                double dy = y - 32.0;
                double dz = z - 32.0;

                // Distance from the tube center-line in the XY plane
                double len_xy = std::sqrt(dx*dx + dy*dy);
                double tube_dist = std::sqrt(std::pow(len_xy - major_R, 2) + dz*dz);

                if (tube_dist <= minor_r) {
                    torus.data[x * ny * nz + y * nz + z] = true;
                }
            }
        }
    }

    // 5. Create Geometry: "The Phantom Peak Cluster"
    // Four large spheres overlapping heavily in a tight square.
    // This creates a Distance Transform peak in the center (32,32,32)
    // which does NOT correspond to any of the original 4 spheres.
    // The algorithm will likely pick the center first, then back-fill the corners.
    VoxelGrid<bool> cluster(nx, ny, nz, v_size);
    
    std::cout << "[1/3] Generating Phantom Peak Cluster..." << std::endl;

    struct BallData { double x, y, z, r; };
    std::vector<BallData> cluster_spheres = {
        {25.0, 25.0, 32.0, 18.0}, // Top-Left
        {25.0, 39.0, 32.0, 18.0}, // Top-Right
        {39.0, 25.0, 32.0, 18.0}, // Bottom-Left
        {39.0, 39.0, 32.0, 18.0}  // Bottom-Right
        
    };

    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int z = 0; z < nz; ++z) {
                bool is_inside = false;
                for (const auto& b : cluster_spheres) {
                    // Check distance squared
                    double d2 = std::pow(x - b.x, 2) + 
                                std::pow(y - b.y, 2) + 
                                std::pow(z - b.z, 2);
                    if (d2 <= b.r * b.r) {
                        is_inside = true;
                        break;
                    }
                }
                if (is_inside) cluster(x, y, z) = true;
            }
        }
    }


    // 7. Create Geometry: "The Sputnik" (Large Core + Satellites)
    // A large central sphere (Radius 20) swallows the inner volume.
    // 6 smaller spheres (Radius 8) protrude slightly from the cardinal directions.
    //
    // WHY THIS BREAKS ITERATION 1:
    // The Global Maxima of the EDT is at (32,32,32) with value ~20.0.
    // The algorithm MUST pick this center first.
    // It cannot "see" the 6 satellite centers until the central volume is masked.
    VoxelGrid<bool> sputnik(nx, ny, nz, v_size);
    
    std::cout << "[1/3] Generating Sputnik..." << std::endl;

    struct SputnikBall { double x, y, z, r; };
    std::vector<SputnikBall> sputnik_geo;

    // The "Sun" (Core)
    sputnik_geo.push_back({32.0, 32.0, 32.0, 20.0}); 

    // The "Satellites" (Protruding bumps)
    // Placed at distance 24 from center. 
    // Overlap = (20 + 8) - 24 = 4 units of overlap.
    double dist = 24.0;
    double sat_r = 8.0;
    
    sputnik_geo.push_back({32.0 + dist, 32.0, 32.0, sat_r}); // +X
    sputnik_geo.push_back({32.0 - dist, 32.0, 32.0, sat_r}); // -X
    sputnik_geo.push_back({32.0, 32.0 + dist, 32.0, sat_r}); // +Y
    sputnik_geo.push_back({32.0, 32.0 - dist, 32.0, sat_r}); // -Y
    sputnik_geo.push_back({32.0, 32.0, 32.0 + dist, sat_r}); // +Z
    sputnik_geo.push_back({32.0, 32.0, 32.0 - dist, sat_r}); // -Z

    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int z = 0; z < nz; ++z) {
                bool is_inside = false;
                for (const auto& b : sputnik_geo) {
                    double d2 = std::pow(x - b.x, 2) + 
                                std::pow(y - b.y, 2) + 
                                std::pow(z - b.z, 2);
                    if (d2 <= b.r * b.r) {
                        is_inside = true;
                        break;
                    }
                }
                if (is_inside) sputnik(x, y, z) = true;
            }
        }
    }

    // ==================== VALIDATION TEST: NON-CUBIC GRID EDT ==================== //
    // 1. Create Non-Cubic Grid (100 x 30 x 30)
    // If layout logic is wrong, the data will wrap incorrectly here.
    int b_nx = 100;
    int b_ny = 30;
    int b_nz = 30;
    VoxelGrid<bool> bar_grid(b_nx, b_ny, b_nz, 1.0);

    // 2. Fill a Solid Bar in the middle
    // Size: 80 x 20 x 20
    // Expected Max Internal Distance: ~10.0 (Center to surface is 10 voxels)
    for(int x=10; x<90; ++x) {
        for(int y=5; y<25; ++y) {
            for(int z=5; z<25; ++z) {
                bar_grid(x,y,z) = true;
            }
        }
    }

    // 8. Create Geometry: "The Tumor" (Merged Gradient)
    // Sphere A: (32,32,32) R=20
    // Sphere B: (44,32,32) R=15
    //
    // Analysis:
    // Distance between centers = 12.
    // Sphere B center is INSIDE Sphere A (since 12 < 20).
    // Therefore, the EDT at (44,32,32) is NOT a local maximum. 
    // The gradient flows uphill toward (32,32,32).
    // The algorithm will see only ONE peak initially.
    // However, Sphere B protrudes: (12 + 15) = 27 vs A's reach of 20.
    // It sticks out by 7 voxels.
    VoxelGrid<bool> tumor(nx, ny, nz, v_size);
    
    std::cout << "[1/3] Generating Tumor..." << std::endl;

    struct TumorBall { double x, y, z, r; };
    std::vector<TumorBall> tumor_geo = {
        {32.0, 32.0, 32.0, 20.0}, // Host
        {44.0, 32.0, 32.0, 15.0}  // Tumor (Center is inside Host)
    };

    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int z = 0; z < nz; ++z) {
                bool is_inside = false;
                for (const auto& b : tumor_geo) {
                    double d2 = std::pow(x - b.x, 2) + 
                                std::pow(y - b.y, 2) + 
                                std::pow(z - b.z, 2);
                    if (d2 <= b.r * b.r) {
                        is_inside = true;
                        break;
                    }
                }
                if (is_inside) tumor(x, y, z) = true;
            }
        }
    }

    // 9. Create Geometry: "Swiss Cheese" (Concave Voids)
    // A solid block (40x40x40) with 3 spherical holes subtracted.
    // This forces the algorithm to pack spheres around the void.
    VoxelGrid<bool> cheese(nx, ny, nz, v_size);
    
    std::cout << "[1/3] Generating Swiss Cheese..." << std::endl;

    // The solid block bounds
    int box_min = 12, box_max = 52; // 40^3 Cube

    // The "Holes" to subtract
    struct Hole { double x, y, z, r; };
    std::vector<Hole> holes = {
        {20.0, 20.0, 52.0, 12.0}, // Corner bite Top-Left-Front
        {45.0, 45.0, 12.0, 15.0}, // Corner bite Bottom-Right-Back
        {32.0, 32.0, 32.0, 10.0}  // Central Void (Internal cavity)
    };

    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int z = 0; z < nz; ++z) {
                // 1. Must be inside the Box
                bool inside_box = (x >= box_min && x < box_max) &&
                                  (y >= box_min && y < box_max) &&
                                  (z >= box_min && z < box_max);
                
                if (!inside_box) continue;

                // 2. Must NOT be inside any Hole
                bool inside_hole = false;
                for (const auto& h : holes) {
                    double d2 = std::pow(x - h.x, 2) + 
                                std::pow(y - h.y, 2) + 
                                std::pow(z - h.z, 2);
                    if (d2 <= h.r * h.r) {
                        inside_hole = true;
                        break;
                    }
                }

                if (!inside_hole) {
                    cheese(x, y, z) = true;
                }
            }
        }
    }

    save_mesh_to_stl(grid_to_mesh(cheese), "swiss_cheese.stl");



    // 10. Create Geometry: "Ice Cream Cone"
    // Geometry Definition:
    // - Cone Tip: (32, 32, 10)
    // - Cone Base: (32, 32, 50) -> Height 40
    // - Base Radius: 10.0 (Aspect Ratio 4:1)
    // - Sphere Cap: Center (32, 32, 50), Radius 12.0 (The "Scoop")
    VoxelGrid<uint8_t> ice_cream(nx, ny, nz, v_size);
    
    std::cout << "[1/3] Generating Ice Cream Cone..." << std::endl;

    double cone_tip_z = 2.0;
    double cone_base_z = 50.0;
    double cone_base_r = 10.0;
    double scoop_r = 12.0; 
    
    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int z = 0; z < nz; ++z) {
                double dx = x - 32.0;
                double dy = y - 32.0;
                
                // --- Cone Logic ---
                bool in_cone = false;
                if (z >= cone_tip_z && z <= cone_base_z) {
                    // Radius increases linearly from tip to base
                    double height_fraction = (z - cone_tip_z) / (cone_base_z - cone_tip_z);
                    double current_r = height_fraction * cone_base_r;
                    
                    if ((dx*dx + dy*dy) <= (current_r * current_r)) {
                        in_cone = true;
                    }
                }

                // --- Sphere Cap Logic ---
                // Centered exactly at the top of the cone (the base)
                double dz_sphere = z - cone_base_z; 
                bool in_sphere = (dx*dx + dy*dy + dz_sphere*dz_sphere) <= (scoop_r * scoop_r);

                if (in_cone || in_sphere) {
                    ice_cream(x, y, z) = true;
                }
            }
        }
    }
    
    save_voxels_to_npy("ice_cream_cone.npy", ice_cream);
    save_mesh_to_stl(grid_to_mesh(ice_cream), "ice_cream_cone.stl");

    //////////////////////////////// Bar Grid EDT Test ///////////////////////////////

    // 3. Run EDT
    // Note: This calls the function in your datatypes header
    auto bar_edt = bar_grid.distance_transform();

    // 4. Analyze Result
    float max_dist = 0.0f;
    for(float v : bar_edt.data) {
        if(v > max_dist) max_dist = v;
    }

    std::cout << "Grid Size: 100 x 30 x 30" << std::endl;
    std::cout << "Bar Size:   80 x 20 x 20 (Thickness 20)" << std::endl;
    std::cout << "Expected Max Distance: ~10.0" << std::endl;
    std::cout << "Computed Max Distance: " << max_dist << std::endl;

    if (max_dist < 5.0) {
        std::cout << "[FAIL] Max distance is too small! Data is being read Transposed/Scrambled." << std::endl;
    } else {
        std::cout << "[PASS] Distance looks valid." << std::endl;
    }

    /////////////////////////////// End of Bar Grid Test /////////////////////////////// 


    std::cout << "[1/3] Synthetic voxel grid created." << std::endl;

    // 2. Run the Reconstruction Algorithm
    // We target 10 spheres or 95% precision
    std::cout << "[2/3] Running reconstruction..." << std::endl;
    SpherePack single_sp = multisphere_from_voxels(
        single_sphere, 
        4,      // min_center_distance_vox
        1, //max iter
        2,      // min_radius_vox
        0.95,   // precision_target
        10,     // max_spheres
        true    // show_progress
    );

    

    SpherePack double_sp = multisphere_from_voxels(
        double_sphere, 
        4,      // min_center_distance_vox
        1, // max iter
        2,      // min_radius_vox
        0.95,   // precision_target
        10,     // max_spheres
        true    // show_progress
    );

    SpherePack rectangle_sp = multisphere_from_voxels(
        rectangle, 
        4,      // min_center_distance_vox
        1, // max iter
        2,      // min_radius_vox
        0.95,   // precision_target
        10,     // max_spheres
        true    // show_progress
    );

    // --- L-SHAPE RECONSTRUCTION ---
    // Note: We perform this with slightly higher sphere count to allow fitting the corner
    SpherePack l_shape_sp = multisphere_from_voxels(
        l_shape, 
        3,      // min_center_distance_vox (Reduced to allow tighter packing in corner)
        1, // max iter
        2,      // min_radius_vox
        0.95,   // precision_target
        15,     // max_spheres
        true    // show_progress
    );

    // --- TORUS RECONSTRUCTION ---
    // A torus requires many spheres to approximate the curve
    SpherePack torus_sp = multisphere_from_voxels(
        torus, 
        3,      // min_center_distance_vox
        1, // max iter
        2,      // min_radius_vox
        0.90,   // precision_target (Slightly lower to prevent infinite small spheres)
        30,     // max_spheres (Increased significantly for the ring)
        true    // show_progress
    );

    // --- CLUSTER RECONSTRUCTION ---
    // We expect the first sphere to be a "phantom" in the middle, 
    // followed by 4 corrective spheres for the corners.
    SpherePack cluster_sp = multisphere_from_voxels(
        cluster, 
        2,      // min_center_distance_vox (Low, because overlap is tight)
        1, // max iter
        2,      // min_radius_vox
        0.95,   // precision_target
        10,     // max_spheres
        true    // show_progress
    );

    SpherePack sputnik_sp = multisphere_from_voxels(
        sputnik, 
        2,      // min_center_distance_vox 
        1, // max iter
        2,      // min_radius_vox
        0.95,   // precision_target
        10,     // max_spheres
        true    // show_progress
    );

    SpherePack tumor_sp = multisphere_from_voxels(
        tumor, 
        2,    // min_center_distance_vox 
        1 , //maxiter
        2,      // min_radius_vox
        0.95,   // precision_target
        10,     // max_spheres
        true    // show_progress
    );

    SpherePack cheese_sp = multisphere_from_voxels(
        cheese, 
        2,      // min_center_distance_vox 
        1, // max iter
        2,      // min_radius_vox (Keep small to fit in tight corners)
        0.92,   // precision_target (Slightly lower as perfect concave fill is hard)
        1050,     // max_spheres (High count required for concave boundaries)
        true    // show_progress
    );

    SpherePack ice_cream_no_boost = multisphere_from_voxels(
        ice_cream,
        15,      // min_center_distance_vox 
        1, // max iter
        2,      // min_radius_vox (Keep small to fit in tight corners)
        0.92,   // precision_target (Slightly lower as perfect concave fill is hard)
        10,     // max_spheres (High count required for concave boundaries)
        true    // show_progress
    );


    SpherePack ice_cream_boost = multisphere_from_voxels(
        ice_cream,
        15,      // min_center_distance_vox 
        10, // max iter
        2,      // min_radius_vox (Keep small to fit in tight corners)
        0.99,   // precision_target (Slightly lower as perfect concave fill is hard)
        10,     // max_spheres (High count required for concave boundaries)
        true    // show_progress
    );
    


    std::cout << "\nReconstruction Complete!" << std::endl;
    std::cout << "--Single Sphere : \n Spheres found: " << single_sp.num_spheres() << std::endl;
    std::cout << "Max radius: " << single_sp.max_radius() << " units" << std::endl;

    
    std::cout << "--Double Sphere : \n Spheres found: " << double_sp.num_spheres() << std::endl;
    std::cout << "Max radius: " << double_sp.max_radius() << " units" << std::endl;

    std::cout << "--Rectangle : \n Spheres found: " << rectangle_sp.num_spheres() << std::endl;
    std::cout << "Max radius: " << rectangle_sp.max_radius() << " units" << std::endl;

    std::cout << "--L-Shape : \n Spheres found: " << l_shape_sp.num_spheres() << std::endl;
    std::cout << "--Torus : \n Spheres found: " << torus_sp.num_spheres() << std::endl;

    std::cout << "--Phantom Peak Cluster : \n Spheres found: " << cluster_sp.num_spheres() << std::endl;

    std::cout << "--Sputnik : \n Spheres found: " << sputnik_sp.num_spheres() << std::endl;

    std::cout << "--Tumor : \n Spheres found: " << tumor_sp.num_spheres() << std::endl;

    std::cout << "--Swiss Cheese : \n Spheres found: " << cheese_sp.num_spheres() << std::endl;

    std::cout << "--Ice Cream Cone (No Boost) : \n Spheres found: " << ice_cream_no_boost.num_spheres() << std::endl;
    std::cout << "--Ice Cream Cone (With Boost) : \n Spheres found: " << ice_cream_boost.num_spheres() << std::endl;



    // 3. Visualization (Optional based on CMake/HAVE_VTK)
#ifdef HAVE_VTK
    std::cout << "[3/3] Opening VTK Visualization..." << std::endl;
    plot_sphere_pack(single_sp);
    plot_sphere_pack(double_sp);
    plot_sphere_pack(rectangle_sp);
#else
    std::cout << "[3/3] VTK not enabled. Skipping visualization." << std::endl;
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

    export_to_csv(ice_cream_no_boost, "reconstructed_ice_cream_no_boost.csv");
    export_to_vtk(ice_cream_no_boost, "reconstructed_ice_cream_no_boost.vtk");

    export_to_csv(ice_cream_boost, "reconstructed_ice_cream_boost.csv");
    export_to_vtk(ice_cream_boost, "reconstructed_ice_cream_boost.vtk");

    std::cout << "To see the result, export to CSV or recompile with VTK fixed." << std::endl;
#endif

    return 0;
}