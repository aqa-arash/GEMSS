#ifndef GEMSS_CONFIG_HPP
#define GEMSS_CONFIG_HPP


#include <Eigen/Core>

namespace GEMSS {


struct MultisphereConfig {
    // --- Required for mesh reconstruction ---
    int div = 200;               // grid resolution
    int padding = 2;             // grid padding
    float minimum_radius_real = 0.0f; // min sphere radius (real units)
    bool confine_mesh = false;   // confine to mesh boundary

    // --- Peak detection and filtering ---
    float min_center_distance_rel = 1.4f; // min center distance (relative) computed as k <= sqrt(min_radius)
    int search_window = 2;       // neighbor distance for peaks
    float radius_offset_vox = 0.0f; // radius correction
    int min_radius_vox = 4;      // min sphere radius (voxels)


    // --- Convergence & limits ---
    float precision_target = 1.0f; // target precision
    int max_spheres = 0;         // max spheres (0 = unlimited)
    int persistence = 2;         // solver persistence

    // --- Utilities & prior state ---
    bool is_2d = false;
    int compute_physics = 2;     // compute physics properties
    bool is_2d = false;   // treat reconstruction as planar (single z-slice + 2D physics)
    float density = 1.0f;        // material density for mass properties 
    bool prune_isolated_spheres = false; // remove isolated spheres
    bool show_progress = true;   // print progress
    Eigen::MatrixX4f initial_sphere_table = Eigen::MatrixX4f(0, 4); // initial spheres

    // --- Fragmentation setting 
    bool conserve_mass = true;
    float parent_representation_offset = 0.87f; // offset for the voxelization of the parent from spheres -> ensures no new spheres place engulfing the old ones (shrink the target by 1 diagonal voxel)
    float fracture_plane_offset = 0.15f; // offset coefficient for fracture plane during splitting ( x * voxel size, applied to the plane normal direction, larger values creates a gap between fragments to ensure clean splitting - model sawed off at the plane without residual voxels bridging the fragments, but causes mass loss)
    
};

} // namespace GEMSS

#endif // GEMSS_CONFIG_HPP
