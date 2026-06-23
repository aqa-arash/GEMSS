/**
 * @file GEMSS_split.hpp
 * @brief SpherePack and voxel grid splitting utilities for GEMSS.
 *
 * Provides functions to split a SpherePack and its voxel grid by a plane, returning separated groups and relabeling the voxel grid.
 *
 * @author Arash Moradian
 * @date 2026-03-31
 */

#ifndef GEMSS_SPLIT_HPP
#define GEMSS_SPLIT_HPP

#include <vector>
#include <Eigen/Core>

#ifdef HAVE_OPENMP
    #include <omp.h>
#endif

#include "GEMSS_datatypes.hpp"
#include "GEMSS_voxel_processing.hpp"
#include "GEMSS_config.hpp"

namespace GEMSS {


/**
 * @brief Computes the projected planar area of the fracture surface between two regions.
 * @param labeled_grid Multi-labeled voxel grid (1: above, 2: below).
 * @param plane_normal Analytical normal vector of the cutting plane.
 * @return Projected surface area of the fracture plane in physical units.
 */
inline float compute_planar_area(
    const VoxelGrid<uint8_t>& labeled_grid,
    const Eigen::Vector3f& plane_normal) 
{
    float total_area = 0.0f;
    float voxel_area = labeled_grid.voxel_size * labeled_grid.voxel_size;
    
    int nx = labeled_grid.nx();
    int ny = labeled_grid.ny();
    int nz = labeled_grid.nz();
    
    Eigen::Vector3f n_plane = plane_normal.normalized();

    // 6 orthogonal face directions for a voxel
    const int dx[] = {1, -1, 0, 0, 0, 0};
    const int dy[] = {0, 0, 1, -1, 0, 0};
    const int dz[] = {0, 0, 0, 0, 1, -1};
    const float nx_dir[] = {1.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    const float ny_dir[] = {0.0f, 0.0f, 1.0f, -1.0f, 0.0f, 0.0f};
    const float nz_dir[] = {0.0f, 0.0f, 0.0f, 0.0f, 1.0f, -1.0f};

    // Parallel reduction to safely accumulate the float value
    #pragma omp parallel for reduction(+:total_area) collapse(3)
    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int z = 0; z < nz; ++z) {
                // Find voxels belonging to Fragment 1 (above plane)
                if (labeled_grid(x, y, z) == 1) { 
                    for (int i = 0; i < 6; ++i) {
                        int nx_idx = x + dx[i];
                        int ny_idx = y + dy[i];
                        int nz_idx = z + dz[i];
                        
                        // Boundary check
                        if (nx_idx >= 0 && nx_idx < nx && ny_idx >= 0 && ny_idx < ny && nz_idx >= 0 && nz_idx < nz) {
                            // If neighbor belongs to Fragment 2, a fracture interface exists
                            if (labeled_grid(nx_idx, ny_idx, nz_idx) == 2) { 
                                Eigen::Vector3f n_face(nx_dir[i], ny_dir[i], nz_dir[i]);
                                
                                // Project the voxel face onto the analytical cutting plane
                                float projection = std::abs(n_face.dot(n_plane));
                                total_area += voxel_area * projection;
                            }
                        }
                    }
                }
            }
        }
    }
    
    return total_area;
}


/**
 * @brief multisphere_from_splitted_voxelGrid
 * @param labeled_grid Multi-labeled voxel grid.
 * @param sphere_lists Vector of SpherePack for each label.
 * @param config MultisphereConfig.
 * @return std::vector<SpherePack> for each region.
 */
inline std::vector<SpherePack> multisphere_from_splitted_voxelGrid(
    const VoxelGrid<uint8_t>& labeled_grid,
    const std::vector<SpherePack>& sphere_lists,
    const MultisphereConfig& config = MultisphereConfig())
{
    int max_label = 0;
    for (size_t i = 0; i < labeled_grid.data.size(); ++i) {
        if (labeled_grid.data[i] > max_label) max_label = labeled_grid.data[i];
    }
    if (max_label == 0) return {};

    std::vector<SpherePack> result(max_label);
    
    #pragma omp parallel for
    for (int label = 1; label <= max_label; ++label) {
        
        // 1. Extract the specific region mask for this fragment
        VoxelGrid<uint8_t> region_grid(labeled_grid.nx(), labeled_grid.ny(), labeled_grid.nz(), labeled_grid.voxel_size, labeled_grid.origin);
        for (int x = 0; x < labeled_grid.nx(); ++x) {
            for (int y = 0; y < labeled_grid.ny(); ++y) {
                for (int z = 0; z < labeled_grid.nz(); ++z) {
                    if (labeled_grid(x, y, z) == label) region_grid(x, y, z) = 1;
                }
            }
        }
        
        // 2. Setup the initial table for this fragment
        MultisphereConfig region_config = config;
        region_config.is_2d = (labeled_grid.nz() == 1);
        if (label - 1 < static_cast<int>(sphere_lists.size()) && sphere_lists[label - 1].centers.rows() > 0) {
            Eigen::MatrixX4f table(sphere_lists[label - 1].centers.rows(), 4);
            // Convert physical coordinates to voxel coordinates for the solver
            table.leftCols(3) = (sphere_lists[label - 1].centers.rowwise() - labeled_grid.origin.transpose()).array() / labeled_grid.voxel_size - 0.5f;
            table.col(3) = sphere_lists[label - 1].radii.array() / labeled_grid.voxel_size;
            region_config.initial_sphere_table = table;
        }
        
        // 3. Reconstruct (This will naturally fill the volume of the discarded boundary spheres)
        result[label - 1] = multisphere_from_voxels(region_grid, region_config);
    }
    
    return result;
}

/**
 * @brief Splits a SpherePack and its voxel grid by a plane defined by a normal and a point.
 *
 * @param sp Input SpherePack.
 * @param normal Plane normal vector (should be normalized).
 * @param point Point on the plane (default: origin).
 * @param config MultisphereConfig for voxelization.
 * @return pair: (List of spherepacks, labeled_voxel_grid)
 */
inline std::pair<std::vector<SpherePack>, VoxelGrid<uint8_t>>
split_sp(const SpherePack& sp, const Eigen::Vector3f& normal, const Eigen::Vector3f& point = Eigen::Vector3f::Zero(), const MultisphereConfig& config = MultisphereConfig()) {
    
    // --- Step 1: Compute bounds and voxelize all spheres ---
    Eigen::Vector3f min_corner = sp.centers.row(0).transpose() - Eigen::Vector3f::Ones() * sp.radii(0);
    Eigen::Vector3f max_corner = sp.centers.row(0).transpose() + Eigen::Vector3f::Ones() * sp.radii(0);
    for (int i = 1; i < sp.centers.rows(); ++i) {
        Eigen::Vector3f c = sp.centers.row(i).transpose();
        float r = sp.radii(i);
        min_corner = min_corner.cwiseMin(c - Eigen::Vector3f::Ones() * r);
        max_corner = max_corner.cwiseMax(c + Eigen::Vector3f::Ones() * r);
    }
    
    float voxel_size = (max_corner - min_corner).maxCoeff() / static_cast<float>(config.div > 0 ? config.div : 100);
    if (voxel_size <= 0) voxel_size = 1.0f;
    Eigen::Vector3f origin = min_corner.array() - config.padding * voxel_size;
    Eigen::Array3i dims = ((max_corner - min_corner).array() / voxel_size + 2 * config.padding).cast<int>();

    if (sp.is_2d) {
        origin.z() = sp.centers(0, 2) - 0.5f * voxel_size;
        dims[2] = 1;
    }
    VoxelGrid<uint8_t> grid(std::max(1, dims[0]), std::max(1, dims[1]),
                            std::max(1, dims[2]), voxel_size, origin);
    
    Eigen::MatrixX4f spheres(sp.centers.rows(), 4);
    spheres.leftCols(3) = (sp.centers.rowwise() - origin.transpose()).array() / voxel_size - 0.5f;
    spheres.col(3) = sp.radii.array() / voxel_size;
    spheres_to_grid(grid, spheres, static_cast<uint8_t>(1));

    // --- Step 2: Compare voxel grid to cutting surface and relabel ---
    VoxelGrid<uint8_t> labeled_grid(grid.nx(), grid.ny(), grid.nz(), grid.voxel_size, grid.origin);
    for (int x = 0; x < grid.nx(); ++x) {
        for (int y = 0; y < grid.ny(); ++y) {
            for (int z = 0; z < grid.nz(); ++z) {
                if (grid(x, y, z) > 0) {
                    // Compute exact physical center of the voxel for accurate splitting
                    Eigen::Vector3f pos = grid.origin + grid.voxel_size * Eigen::Vector3f(x + 0.5f, y + 0.5f, z + 0.5f);
                    float d = normal.dot(pos - point);
                    // >= ensures voxels exactly on the boundary aren't lost
                    if (d >= 0) labeled_grid(x, y, z) = 1; 
                    else labeled_grid(x, y, z) = 2;
                }
            }
        }
    }

    // --- Step 3: Compute exact sphere contact and group non-intersecting spheres ---
    std::vector<std::vector<int>> group_indices(2); // 0: above (label 1), 1: below (label 2)
    for (int i = 0; i < sp.centers.rows(); ++i) {
        float d = normal.dot(sp.centers.row(i).transpose() - point);
        float r = sp.radii(i);
        
        if (d > r) {
            group_indices[0].push_back(i); // Strictly above plane
        } else if (d < -r) {
            group_indices[1].push_back(i); // Strictly below plane
        }
        // Intersecting spheres (|d| <= r) are entirely ignored and deleted
    }
    
    std::vector<SpherePack> sphere_groups;
    for (int g = 0; g < 2; ++g) {
        Eigen::MatrixX3f centers(group_indices[g].size(), 3);
        Eigen::VectorXf radii(group_indices[g].size());
        for (size_t i = 0; i < group_indices[g].size(); ++i) {
            centers.row(i) = sp.centers.row(group_indices[g][i]);
            radii(i) = sp.radii(group_indices[g][i]);
        }
        sphere_groups.emplace_back(centers, radii);
    }

    // --- Step 4: Pass pairs to build function ---
    std::vector<SpherePack> reconstructed = multisphere_from_splitted_voxelGrid(labeled_grid, sphere_groups, config);

     // --- Step 5: Strict Mass Conservation ---
    if (config.conserve_mass && sp.mass > 0.0f) 
    {
        float total_new_mass = 0.0f;
        // Calculate the uncorrected total mass of all generated fragments
        for (const auto& frag : reconstructed) {
            total_new_mass += frag.mass;
        }
        
        if (total_new_mass > 0.0f) 
        {
            // Calculate the scaling coefficient required to strictly match the parent
            float mass_scale = sp.mass / total_new_mass;
            
            // Apportion the parent's mass based on the fragment's volume percentage
            for (auto& frag : reconstructed) {
                frag.mass *= mass_scale;
                frag.density *= mass_scale;             // Linearly scales mass
                frag.inertia_tensor *= mass_scale;      // Linearly scales moments
                frag.principal_moments *= mass_scale;   // Linearly scales moments
                
                // Note: Center of mass and principal axes (eigenvectors) remain mathematically
                // unchanged when density is scaled uniformly, so we do not touch them.
            }
        }
    }
    
    return std::make_pair(reconstructed, labeled_grid);
}


/**
 * @brief Fused loop: Splits a SpherePack and concurrently computes the projected fracture area.
 *
 * @param sp Input SpherePack.
 * @param normal Plane normal vector (should be normalized).
 * @param point Point on the plane (default: origin).
 * @param config MultisphereConfig for voxelization.
 * @return Tuple: (vector of reconstructed SpherePacks, labeled_voxel_grid, projected_fracture_area)
 */
inline std::tuple<std::vector<SpherePack>, VoxelGrid<uint8_t>, float>
split_and_compute_surface_sp(const SpherePack& sp, const Eigen::Vector3f& normal, const Eigen::Vector3f& point = Eigen::Vector3f::Zero(), const MultisphereConfig& config = MultisphereConfig()) {
    
    // --- Steo 0 -> set the density of the new particles equal to their parent particle 
    MultisphereConfig local_config = config;
    local_config.density = sp.density;


    // --- Step 1: Compute bounds and voxelize all spheres ---
    Eigen::Vector3f min_corner = sp.centers.row(0).transpose() - Eigen::Vector3f::Ones() * sp.radii(0);
    Eigen::Vector3f max_corner = sp.centers.row(0).transpose() + Eigen::Vector3f::Ones() * sp.radii(0);
    for (int i = 1; i < sp.centers.rows(); ++i) {
        Eigen::Vector3f c = sp.centers.row(i).transpose();
        float r = sp.radii(i);
        min_corner = min_corner.cwiseMin(c - Eigen::Vector3f::Ones() * r);
        max_corner = max_corner.cwiseMax(c + Eigen::Vector3f::Ones() * r);
    }

    float voxel_size = (max_corner - min_corner).maxCoeff() / static_cast<float>(local_config.div);
    if (voxel_size <= 0) voxel_size = 1.0f;
    Eigen::Vector3f origin = min_corner.array() - local_config.padding * voxel_size;
    Eigen::Array3i dims = ((max_corner - min_corner).array() / voxel_size + 2 * local_config.padding).cast<int>();

    if (sp.is_2d) {
        origin.z() = sp.centers(0, 2) - 0.5f * voxel_size;
        dims[2] = 1;
    }
    VoxelGrid<uint8_t> grid(std::max(1, dims[0]), std::max(1, dims[1]),
                            std::max(1, dims[2]), voxel_size, origin);
    
    Eigen::MatrixX4f spheres(sp.centers.rows(), 4);
    spheres.leftCols(3) = (sp.centers.rowwise() - origin.transpose()).array() / voxel_size - 0.5f;
    spheres.col(3) = sp.radii.array() / voxel_size - local_config.parent_representation_offset;
    spheres_to_grid(grid, spheres, static_cast<uint8_t>(1));

    // --- Step 2: Fused Labeling and Surface Area Projection ---
    VoxelGrid<uint8_t> labeled_grid(grid.nx(), grid.ny(), grid.nz(), grid.voxel_size, grid.origin);
    
    float total_area = 0.0f;
    float voxel_area = grid.voxel_size * grid.voxel_size;
    Eigen::Vector3f n_plane = normal.normalized();

    // 6 orthogonal face directions for a voxel
    const int dx[] = {1, -1, 0, 0, 0, 0};
    const int dy[] = {0, 0, 1, -1, 0, 0};
    const int dz[] = {0, 0, 0, 0, 1, -1};
    const float nx_dir[] = {1.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    const float ny_dir[] = {0.0f, 0.0f, 1.0f, -1.0f, 0.0f, 0.0f};
    const float nz_dir[] = {0.0f, 0.0f, 0.0f, 0.0f, 1.0f, -1.0f};

    // Parallelize with reduction for the float area accumulator
    #pragma omp parallel for reduction(+:total_area) collapse(3)
    for (int x = 0; x < grid.nx(); ++x) {
        for (int y = 0; y < grid.ny(); ++y) {
            for (int z = 0; z < grid.nz(); ++z) {
                if (grid(x, y, z) > 0) { // Only process parent solid voxels
                    
                    Eigen::Vector3f pos = grid.origin + grid.voxel_size * Eigen::Vector3f(x + 0.5f, y + 0.5f, z + 0.5f);
                    float d = n_plane.dot(pos - point);
                    
                    if (d >= local_config.fracture_plane_offset * grid.voxel_size) {
                        labeled_grid(x, y, z) = 1; 
                        
                        // FUSED: Evaluate 6 neighbors mathematically
                        for (int i = 0; i < 6; ++i) {
                            int nx_idx = x + dx[i];
                            int ny_idx = y + dy[i];
                            int nz_idx = z + dz[i];
                            
                            // Check grid bounds and verify neighbor is solid (part of the parent grid)
                            if (nx_idx >= 0 && nx_idx < grid.nx() && 
                                ny_idx >= 0 && ny_idx < grid.ny() && 
                                nz_idx >= 0 && nz_idx < grid.nz() && 
                                grid(nx_idx, ny_idx, nz_idx) > 0) 
                            {
                                Eigen::Vector3f n_face(nx_dir[i], ny_dir[i], nz_dir[i]);
                                // Exact mathematical distance of the neighbor to the plane
                                float d_neighbor = d + grid.voxel_size * n_plane.dot(n_face);
                                
                                // If the solid neighbor lies across the cut plane, a new surface is born
                                if (d_neighbor < 0) {
                                    total_area += voxel_area * std::abs(n_face.dot(n_plane));
                                }
                            }
                        }
                    } else if (d <= - local_config.fracture_plane_offset * grid.voxel_size){
                        labeled_grid(x, y, z) = 2;
                    } else {
                        labeled_grid(x,y,z) = 0;
                    }
                    
                }
            }
        }
    }

    // --- Step 3: Compute exact sphere contact and group non-intersecting spheres ---
    std::vector<std::vector<int>> group_indices(2);
    for (int i = 0; i < sp.centers.rows(); ++i) {
        float d = normal.dot(sp.centers.row(i).transpose() - point);
        float r = sp.radii(i);
        
        if (d > r) {
            group_indices[0].push_back(i); 
        } else if (d < -r) {
            group_indices[1].push_back(i); 
        }
    }
    
    std::vector<SpherePack> sphere_groups;
    for (int g = 0; g < 2; ++g) {
        Eigen::MatrixX3f centers(group_indices[g].size(), 3);
        Eigen::VectorXf radii(group_indices[g].size());
        for (size_t i = 0; i < group_indices[g].size(); ++i) {
            centers.row(i) = sp.centers.row(group_indices[g][i]);
            radii(i) = sp.radii(group_indices[g][i]);
        }
        sphere_groups.emplace_back(centers, radii);
    }

    // --- Step 4: Pass pairs to build function ---
    std::vector<SpherePack> reconstructed = multisphere_from_splitted_voxelGrid(labeled_grid, sphere_groups, local_config);

    // --- Step 5: Strict Mass Conservation ---
    if (local_config.conserve_mass && sp.mass > 0.0f) 
    {
        float total_new_mass = 0.0f;
        // Calculate the uncorrected total mass of all generated fragments
        for (const auto& frag : reconstructed) {
            total_new_mass += frag.mass;
        }
        
        if (total_new_mass > 0.0f) 
        {
            // Calculate the scaling coefficient required to strictly match the parent
            float mass_scale = sp.mass / total_new_mass;
            
            // Apportion the parent's mass based on the fragment's volume percentage
            for (auto& frag : reconstructed) {
                frag.mass *= mass_scale;
                frag.density *= mass_scale;             // Linearly scales mass
                frag.inertia_tensor *= mass_scale;      // Linearly scales moments
                frag.principal_moments *= mass_scale;   // Linearly scales moments
                
                // Note: Center of mass and principal axes (eigenvectors) remain mathematically
                // unchanged when density is scaled uniformly, so we do not touch them.
            }
        }
    }
    
    return std::make_tuple(reconstructed, labeled_grid, total_area);
}




} // namespace GEMSS

#endif // GEMSS_SPLIT_HPP