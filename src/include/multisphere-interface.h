#ifndef MULTISPHERE_INTERFACE_H
#define MULTISPHERE_INTERFACE_H

/**
 * @file multisphere-interface.h
 * @brief Umbrella header for the multisphere-cpp library.
 *
 * Includes all main data structures and user-facing API functions for mesh voxelization,
 * multisphere reconstruction, I/O, and voxel processing. Include this file to access the
 * primary functionality of the library.
 *
 * @author Arash Moradian
 * @date 2026-03-10
 */




#include "multisphere_datatypes.hpp"
#include "multisphere_io.hpp"
#include "multisphere_mesh_handler.hpp"
#include "multisphere_reconstruction_helpers.hpp"
#include "multisphere_reconstruction.hpp"
#include "multisphere_voxel_processing.hpp"
#include "multisphere_config.hpp"

namespace MSS {

// --- Main Data Structures ---
// See multisphere_datatypes.hpp for details
// (Types and templates are now in MSS namespace)

// --- Main API Functions ---

/**
 * @brief Load a mesh from a binary STL file.
 * @param path Path to STL file.
 * @return FastMesh structure.
 */
inline FastMesh load_mesh_fast(const std::string& path);

/**
 * @brief Save a FastMesh to a binary STL file.
 * @param mesh FastMesh to save.
 * @param output_path Output STL file path.
 */
inline void save_mesh_to_stl(const FastMesh& mesh, const std::string& output_path);

/**
 * @brief Voxelize a mesh using the robust winding number method.
 * @param mesh Input FastMesh.
 * @param div Voxel grid division (resolution).
 * @param padding Grid padding.
 * @return VoxelGrid<bool> representing mesh occupancy.
 */
inline VoxelGrid<uint8_t> mesh_to_binary_grid(const FastMesh& mesh, int div, int padding /* = 2 */);

/**
 * @brief Multisphere reconstruction from a voxel grid.
 * @tparam T VoxelGrid data type.
 * @param input_grid Input voxel grid.
 * @param config MultisphereConfig struct containing all configuration parameters:
 * @param min_center_distance_vox Minimum center distance (voxels).
 * @param min_radius_vox Minimum radius (voxels).
 * @param precision_target Target precision.
 * @param max_spheres Maximum number of spheres.
 * @param show_progress Show progress output.
 * @param sphere_table_in Optional initial sphere table.
 * @param compute_physics Whether to compute physical properties of the multisphere union. 0 = false, 1 = compute based on reconstruction
 * @return SpherePack reconstruction result.
 */
template <typename T>
SpherePack multisphere_from_voxels(
	const VoxelGrid<T>& input_grid,
	const MultisphereConfig& config //= MultisphereConfig()
);

/**
 * @brief Multisphere reconstruction directly from a triangle mesh.
 * @param mesh Input FastMesh.
 * @param config MultisphereConfig struct containing all configuration parameters:
 * @param div Voxel grid division (resolution).
 * @param padding Grid padding.
 * @param min_center_distance_vox Minimum center distance (voxels).
 * @param min_radius_vox Minimum radius (voxels).
 * @param precision Target precision.
 * @param max_spheres Maximum number of spheres.
 * @param show_progress Show progress output.
 * @param confine_mesh Confine spheres to mesh boundary.
 * @param sphere_table Optional initial sphere table.
 * @param compute_physics Whether to compute physical properties of the multisphere union. 0 = false, 1 = compute based on reconstruction, 2 = compute based on original mesh (if available)
 * @return SpherePack reconstruction result.
 */
SpherePack multisphere_from_mesh(
	const FastMesh& mesh,
	const MultisphereConfig& config //= MultisphereConfig()
);

/**
 * @brief Export a SpherePack to CSV.
 * @param sp SpherePack to export.
 * @param path Output CSV file path.
 */
inline void export_to_csv(const SpherePack& sp, const std::string& path);

/**
 * @brief Export a SpherePack to legacy VTK format.
 * @param sp SpherePack to export.
 * @param path Output VTK file path.
 */
inline void export_to_vtk(const SpherePack& sp, const std::string& path);

/**
 * @brief Export a VoxelGrid distance transform to a legacy VTK file.
 * @tparam T VoxelGrid data type.
 * @param grid VoxelGrid to export.
 * @param path Output VTK file path.
 */
template <typename T>
inline void export_voxel_grid_to_vtk(const VoxelGrid<T>& grid, const std::string& path);

// . Mesh handling utilities .

/** 
 * @brief Compute the minimum axis-aligned bounding box (AABB) dimension of a mesh.
 * @details Used for converting from voxel units to physical units and for setting default parameters based on mesh scale.
 * @param mesh Input FastMesh.
 * @return Minimum AABB dimension (float).
 */
inline float get_min_AABB(const FastMesh & mesh);

/**
 * @brief Filters the sphere table to keep only the largest connected network.
 * @details Used to ensure a continuous representation of the multisphere union.
 * Uses an OpenMP-parallelized Dense Grid Spatial Hash with dynamic search extents.
 * @param sphere_table Input sphere table (Nx4 matrix).
 * @return Filtered sphere table containing only the primary connected component.
 */
inline Eigen::MatrixX4f filter_largest_sphere_network(const Eigen::MatrixX4f& sphere_table);

} // namespace MSS



#endif // MULTISPHERE_INTERFACE_H