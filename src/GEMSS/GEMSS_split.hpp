
#ifndef GEMSS_SPLIT_HPP
#define GEMSS_SPLIT_HPP

/**
 * @file GEMSS_split.hpp
 * @brief SpherePack and voxel grid splitting utilities for GEMSS.
 *
 * Provides functions to split a SpherePack and its voxel grid by a plane, returning separated groups and relabeling the voxel grid.
 *
 * @author Arash Moradian
 * @date 2026-03-31
 */

#include <vector>
#include <Eigen/Core>
#include "GEMSS_datatypes.hpp"
#include "GEMSS_voxel_processing.hpp"
#include "GEMSS_config.hpp"

namespace GEMSS {

/**
 * @brief Splits a SpherePack and its voxel grid by a plane defined by a normal and a point.
 *
 * @param sp Input SpherePack.
 * @param normal Plane normal vector (should be normalized).
 * @param point Point on the plane (default: origin).
 * @param config MultisphereConfig for voxelization.
 * @return Tuple: (spheres_above, spheres_below, labeled_voxel_grid)
 */

inline std::pair<std::vector<SpherePack>, VoxelGrid<uint8_t>>
split_sp(const SpherePack& sp, const Eigen::Vector3f& normal, const Eigen::Vector3f& point = Eigen::Vector3f::Zero(), const MultisphereConfig& config = MultisphereConfig()) {
    // 1. Compute bounding box considering radii
    Eigen::Vector3f min_corner = sp.centers.row(0).transpose() - Eigen::Vector3f::Ones() * sp.radii(0);
    Eigen::Vector3f max_corner = sp.centers.row(0).transpose() + Eigen::Vector3f::Ones() * sp.radii(0);
    for (int i = 1; i < sp.centers.rows(); ++i) {
        Eigen::Vector3f c = sp.centers.row(i).transpose();
        float r = sp.radii(i);
        min_corner = min_corner.cwiseMin(c - Eigen::Vector3f::Ones() * r);
        max_corner = max_corner.cwiseMax(c + Eigen::Vector3f::Ones() * r);
    }
    float voxel_size = (max_corner - min_corner).minCoeff() / static_cast<float>(config.div > 0 ? config.div : 100);
    if (voxel_size <= 0) voxel_size = 1.0f;
    Eigen::Vector3f origin = min_corner.array() - config.padding * voxel_size;
    Eigen::Array3i dims = ((max_corner - min_corner) / voxel_size).cast<int>() + 2 * config.padding;
    VoxelGrid<uint8_t> grid(dims[0], dims[1], dims[2], voxel_size, origin);
    Eigen::MatrixX4f spheres(sp.centers.rows(), 4);
    spheres.leftCols(3) = sp.centers;
    spheres.col(3) = sp.radii;
    spheres_to_grid(grid, spheres, 1.0f);

    // 2. Split spheres by plane (generalized for more fragments if needed)
    std::vector<std::vector<int>> group_indices(2); // 0: above, 1: below
    for (int i = 0; i < sp.centers.rows(); ++i) {
        float d = normal.dot(sp.centers.row(i).transpose() - point);
        if (d > 0) group_indices[0].push_back(i); // above
        else group_indices[1].push_back(i); // below
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

    // 3. Split voxel grid by plane (label: 1=above, 2=below, 0=on plane)
    VoxelGrid<uint8_t> labeled_grid(grid.nx(), grid.ny(), grid.nz(), grid.voxel_size, grid.origin);
    for (int x = 0; x < grid.nx(); ++x) {
        for (int y = 0; y < grid.ny(); ++y) {
            for (int z = 0; z < grid.nz(); ++z) {
                if (grid(x, y, z) > 0) {
                    Eigen::Vector3f pos = grid.origin + grid.voxel_size * Eigen::Vector3f(x, y, z);
                    float d = normal.dot(pos - point);
                    if (d > 0) labeled_grid(x, y, z) = 1;
                    else if (d < 0) labeled_grid(x, y, z) = 2;
                    else labeled_grid(x, y, z) = 0;
                }
            }
        }
    }

    // 4. Call multisphere_from_splitted_voxelGrid with the new API
    std::vector<SpherePack> reconstructed = multisphere_from_splitted_voxelGrid(labeled_grid, sphere_groups, config);
    return std::make_pair(reconstructed, labeled_grid);
}

/**
 * @brief Placeholder for multisphere_from_splitted_voxelGrid. To be implemented by user.
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
	// Find max label (number of regions)
	int max_label = 0;
	for (size_t i = 0; i < labeled_grid.data.size(); ++i) {
		if (labeled_grid.data[i] > max_label) max_label = labeled_grid.data[i];
	}
	if (max_label == 0) return {};

	std::vector<SpherePack> result(max_label);
	#ifdef HAVE_OPENMP
	#pragma omp parallel for
	#endif
	for (int label = 1; label <= max_label; ++label) {
		// Extract region for this label
		VoxelGrid<uint8_t> region_grid(labeled_grid.nx(), labeled_grid.ny(), labeled_grid.nz(), labeled_grid.voxel_size, labeled_grid.origin);
		for (int x = 0; x < labeled_grid.nx(); ++x) {
			for (int y = 0; y < labeled_grid.ny(); ++y) {
				for (int z = 0; z < labeled_grid.nz(); ++z) {
					if (labeled_grid(x, y, z) == label) region_grid(x, y, z) = 1;
				}
			}
		}
		MultisphereConfig region_config = config;
		if (label-1 < static_cast<int>(sphere_lists.size()) && sphere_lists[label-1].centers.rows() > 0 && sphere_lists[label-1].radii.size() > 0) {
			Eigen::MatrixX4f table(sphere_lists[label-1].centers.rows(), 4);
			table.leftCols(3) = sphere_lists[label-1].centers;
			table.col(3) = sphere_lists[label-1].radii;
			region_config.initial_sphere_table = table;
		}
		result[label-1] = multisphere_from_voxels(region_grid, region_config);
	}
	return result;
}

} // namespace GEMSS

#endif // GEMSS_SPLIT_HPP
