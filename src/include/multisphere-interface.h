
/**
 * @file multisphere-cpp.h
 * @brief Umbrella header for the multisphere-cpp library.
 *
 * Includes all main data structures and user-facing API functions for mesh voxelization,
 * multisphere reconstruction, I/O, and voxel processing. Include this file to access the
 * primary functionality of the library.
 *
 * @author Arash Moradian
 * @date 2026-03-10
 */


#pragma once

/**
 * @file multisphere-cpp.h
 * @brief Umbrella header for the multisphere-cpp library (namespace MSS).
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
inline VoxelGrid<bool> mesh_to_binary_grid(const FastMesh& mesh, int div, int padding /* = 2 */);

/**
 * @brief Multisphere reconstruction from a voxel grid.
 * @tparam T VoxelGrid data type.
 * @param input_grid Input voxel grid.
 * @param min_center_distance_vox Minimum center distance (voxels).
 * @param max_iter Maximum iterations.
 * @param min_radius_vox Minimum radius (voxels).
 * @param precision_target Target precision.
 * @param max_spheres Maximum number of spheres.
 * @param show_progress Show progress output.
 * @param sphere_table_in Optional initial sphere table.
 * @return SpherePack reconstruction result.
 */
template <typename T>
SpherePack multisphere_from_voxels(
	const VoxelGrid<T>& input_grid,
	int min_center_distance_vox, // = 2
	int max_iter, // = 10
	std::optional<int> min_radius_vox, // = std::nullopt
	std::optional<float> precision_target, // = std::nullopt
	std::optional<int> max_spheres, // = std::nullopt
	bool show_progress, // = true
	std::optional<Eigen::MatrixX4f> sphere_table_in // = std::nullopt
);

/**
 * @brief Multisphere reconstruction directly from a triangle mesh.
 * @param mesh Input FastMesh.
 * @param div Voxel grid division (resolution).
 * @param padding Grid padding.
 * @param min_center_distance_vox Minimum center distance (voxels).
 * @param max_iter Maximum iterations.
 * @param min_radius_vox Minimum radius (voxels).
 * @param precision Target precision.
 * @param max_spheres Maximum number of spheres.
 * @param show_progress Show progress output.
 * @param confine_mesh Confine spheres to mesh boundary.
 * @param sphere_table Optional initial sphere table.
 * @return SpherePack reconstruction result.
 */
SpherePack multisphere_from_mesh(
	const FastMesh& mesh,
	int div, // = 100
	int padding, // = 2
	int min_center_distance_vox, // = 4
	int max_iter, // = 1
	std::optional<int> min_radius_vox, // = std::nullopt
	std::optional<float> precision, // = std::nullopt
	std::optional<int> max_spheres, // = std::nullopt
	bool show_progress, // = true
	bool confine_mesh, // = false
	std::optional<Eigen::MatrixX4f> sphere_table // = std::nullopt
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

/**
 * @brief Load a VoxelGrid<bool> from a .npy file.
 * @param path Path to .npy file.
 * @param voxel_size Physical voxel size.
 * @return VoxelGrid<bool> loaded from file.
 */
inline VoxelGrid<bool> load_voxels_from_npy(const std::string& path, double voxel_size);

/**
 * @brief Save a VoxelGrid to a .npy file.
 * @tparam T VoxelGrid data type.
 * @param path Output .npy file path.
 * @param grid VoxelGrid to save.
 */
template<typename T>
inline void save_voxels_to_npy(const std::string& path, const VoxelGrid<T>& grid);

// ... Add more user-facing API declarations as needed ...

/** 
 * @brief Compute the minimum axis-aligned bounding box (AABB) dimension of a mesh.
 * @param mesh Input FastMesh.
 * @return Minimum AABB dimension (float).
 */
inline float get_min_AABB(const FastMesh & mesh);

} // namespace MSS



