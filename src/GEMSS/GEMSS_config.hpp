#ifndef GEMSS_CONFIG_HPP
#define GEMSS_CONFIG_HPP


#include <Eigen/Core>

namespace GEMSS {

/**
 * @brief Configuration parameters for multisphere reconstruction.
 */
struct MultisphereConfig {
    // --- Required parameters for reconstruction from mesh ---
    int div = 100;               ///< Voxel grid division (resolution)
    int padding = 2;             ///< Voxel grid padding
    float minimum_radius_real = 0.0f; ///< Minimum sphere radius in real units (e.g., millimeters). If not set, Voxel units will be used and the min_radius_vox parameter will be required.
    bool confine_mesh = false;   ///< Confine spheres strictly to the mesh boundary SDF

    // --- Peak Detection & Filtering ---
    float min_center_distance_rel = 0.5f;              ///< Minimum center distance relative to radius (e.g., 0.5 means centers must be at least half a radius apart)

    // --- Convergence & Limits ---
    int min_radius_vox = 2; ///< Minimum sphere radius in voxels default is 2 voxels, which allows for capturing small features while avoiding excessive numbers of tiny spheres due to discretization error. Adjust based on the expected feature size in your input data and desired level of detail.
    float precision_target = 1.0f; ///< Target voxel overlap precision [0, 1]. please note 1 is impossible to reach due to discretization, set to 0.99 or lower for practical purposes. Adjust based on the desired balance between reconstruction quality and runtime.
    int max_spheres = 0;        ///< Maximum allowed spheres in the reconstruction (0 for unlimited)

    // --- Utilities & Prior State ---
    int compute_physics = 2; ///< Compute volume, CoM, and inertia tensor 0 = false, 1 = Compute based on reconstruction, 2 = compute based on original mesh (if available)
    bool prune_isolated_spheres = false; ///< Remove spheres that are not touching the biggest network of spheres
    bool show_progress = true;    ///< Print console progress
    Eigen::MatrixX4f initial_sphere_table = Eigen::MatrixX4f(0, 4); ///< Prior solver state, can be in voxel units (if passed to multisphere_from_voxels) or physical units (if passed to multisphere_from_mesh). This is a matrix of shape (N, 4) where each row is (center_x, center_y, center_z, radius). If provided, the solver will start with this initial configuration instead of an empty table. Adjust based on whether you have a good initial guess for the sphere configuration that can speed up convergence.


    // --- Magic numbers (Don't touch unless you know what you're doing) ---
    float radius_offset_vox = 0.87f; ///< Correction factor added to radius during peak detection to better fit the original shape. Default is approximately the half-diagonal of a unit voxel sqrt(3)/2, which empirically improves reconstruction quality. Adjust with caution. 
    int search_window = 2; ///< Minimum distance of neighbors in voxels to consider during peak detection. This is a hard threshold that can be used to smooth out small fluctuation in the distance field and reduce the number of spheres, but setting it too high may cause missed detections in thin structures. Default is 2 voxels, to avoid discretization noise while still allowing for close peaks in thin structures. Adjust based on the voxel resolution and expected feature size in the input data.
    int persistence = 2; ///< Maximum number of iterations for the solver to increase the weight before giving up. Used only when precision_target or max_sphere is not reached. Can cause excessive runtime, set to 1 to turn off. Default is 2, which allows for increasing the weight if no peak is found for 1 iteration, but can be increased for more challenging reconstructions at the cost of runtime. Adjust based on the complexity of your input data and desired precision.

};

} // namespace MSS

#endif // GEMSS_CONFIG_HPP
