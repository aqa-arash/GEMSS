#ifndef MULTISPHERE_CONFIG_HPP
#define MULTISPHERE_CONFIG_HPP

#include <optional>
#include <Eigen/Dense>

namespace MSS {

/**
 * @brief Configuration parameters for multisphere reconstruction.
 */
struct MultisphereConfig {
    // --- Grid Generation (Used only by multisphere_from_mesh) ---
    int div = 100;               ///< Voxel grid division (resolution)
    int padding = 2;             ///< Voxel grid padding
    bool confine_mesh = false;   ///< Confine spheres strictly to the mesh boundary SDF

    // --- Peak Detection & Filtering ---
    int min_center_distance_vox = 2;                  ///< Minimum center distance in voxels
    std::optional<int> min_radius_vox = std::nullopt; ///< Minimum sphere radius in voxels

    // --- Convergence & Limits ---
    std::optional<float> precision_target = std::nullopt; ///< Target voxel overlap precision [0, 1]
    std::optional<int> max_spheres = std::nullopt;        ///< Maximum allowed spheres

    // --- Utilities & Prior State ---
    int compute_physics = 0; ///< Compute volume, CoM, and inertia tensor 0 = false, 1 = Compute based on reconstruction, 2 = compute based on original mesh (if available)
    bool prune_isolated_spheres = false; ///< Remove spheres that are not touching the biggest network of spheres
    bool show_progress = true;    ///< Print console progress
    std::optional<Eigen::MatrixX4f> initial_sphere_table = std::nullopt; ///< Prior solver state
};

} // namespace MSS

#endif // MULTISPHERE_CONFIG_HPP