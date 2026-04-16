# CONFIG_GUIDE.md: MultisphereConfig Structure

This document explains the `MultisphereConfig` structure found in `src/GEMSS/GEMSS_config.hpp`. This structure contains configuration parameters for multisphere reconstruction in the GEMSS library.

## Structure Definition

```cpp
struct MultisphereConfig {
    // --- Required for mesh reconstruction ---
    int div = 100;               // grid resolution
    int padding = 2;             // grid padding
    float minimum_radius_real = 0.0f; // min sphere radius (real units)
    bool confine_mesh = false;   // confine to mesh boundary

    // --- Peak detection and filtering ---
    float min_center_distance_rel = 0.5f; // min center distance (relative)
    int search_window = 2;       // neighbor distance for peaks
    float radius_offset_vox = 0.5f; // radius correction


    // --- Convergence & limits ---
    int min_radius_vox = 2;      // min sphere radius (voxels)
    float precision_target = 1.0f; // target precision
    int max_spheres = 0;         // max spheres (0 = unlimited)
    int persistence = 2;         // solver persistence

    // --- Utilities & prior state ---
    int compute_physics = 2;     // compute physics properties
    bool prune_isolated_spheres = false; // remove isolated spheres
    bool show_progress = true;   // print progress
    Eigen::MatrixX4f initial_sphere_table = Eigen::MatrixX4f(0, 4); // initial spheres
};
```

## Field Descriptions

### Required parameters for reconstruction from mesh
* **div**: Voxel grid division. This determines the resolution by computing a cubic cell size, calculated as `min(AABB) / div` (where AABB is the axis-aligned bounding box). A rectangular grid is then generated based on this cell size.
* **padding**: Padding added to each side of the bounding box grid, measured in voxels.
* **minimum_radius_real**: Minimum sphere radius in real physical units (e.g., millimeters). If set to `0.0f`, physical units are ignored, and the solver falls back to using the `min_radius_vox` parameter.
* **confine_mesh**: Confine spheres strictly to the mesh boundary's Signed Distance Field (SDF).

### Peak Detection & Filtering
  * **min_center_distance_rel**: Minimum center distance relative to the square root of the local radius (acting as coefficient $k$). For example, a value of `0.5` means adjacent centers must be at least $0.5 \times \sqrt{R}$ apart. 
  * **Purpose**: Scales the spacing to maintain a uniform artificial roughness across the shape, preventing over-resolution on flat surfaces while preserving high-curvature details.
  * **Warning (Detachment Limit)**: If the minimum allowable radius in the domain falls below $k^2 / 4$, the required spacing exceeds the diameter of those small spheres, causing them to detach with no physical overlap.
  * **Warning (Choking Limit)**: If the minimum allowable radius falls below $k^2$, the required spacing exceeds the radius itself ($d > R$). In converging geometries (e.g., conical tips), the oversized exclusion zone chokes the domain, and new spheres may fail to generate in these regions.* **search_window**: The number of cells looked at from each direction to assess whether the current cell is a local maximum. The default is `2` to avoid voxelization error. It can be set to higher values to smooth the surface of the geometry. Empirically, a value equal to 4% of the voxelization resolution (`div`) has shown to perform well for complex geometries like sand particle scans.
* **radius_offset_vox**: Correction factor added to the radius during peak detection to better fit the original shape. Default is approximately half a voxel (`0.5f`), which empirically improves reconstruction quality. Adjust with caution.

### Convergence & Limits
* **min_radius_vox**: Minimum sphere radius in voxels. Consider setting it to `2-4%` of the voxelization resolution (`div`), to capture small features while avoiding an excessive number of tiny spheres due to discretization error. Default is `2` voxels for `div = 100`. Will be overwritten if `minimum_radius_real` is set and input is an STL mesh.
* **precision_target**: Target voxel overlap precision `[0, 1]`. **Note on Defaults:** By default, `precision_target` is set to `1.0f` and `max_spheres` to `0`. Because true 1.0 precision is impossible due to discretization, this combination tells the algorithm to run in **exhaustion mode**—meaning it will continue adding spheres until it literally cannot fit any more. If you want faster runtimes, lower this to a practical target (e.g., `0.95` or `0.99`).
* **max_spheres**: Maximum allowed spheres in the reconstruction (`0` = unlimited / run to exhaustion).
* **persistence**: Maximum number of consecutive iterations for the solver to increase the weight without finding new peaks before giving up. Used only when `precision_target` or `max_spheres` is not yet reached. Default is `2`, allowing the solver to increase weight if no peak is found for 1 iteration. Can be increased for challenging reconstructions, but setting it higher will increase runtime. Set to `1` to disable.

### Utilities & Prior State
* **compute_physics**: Determines how volume, center of mass, and inertia tensors are calculated. 
    * `0` = Do not compute.
    * `1` = Compute properties based on the *generated multisphere representation*.
    * `2` = Compute properties based on the *target voxel grid*. (Use this default to create a multisphere approximation while strictly maintaining the physical behavior of the original mesh/voxel grid).
* **prune_isolated_spheres**: If `true`, removes spheres that are not physically touching the largest contiguous network of spheres. Useful when complex geometries are modeled using a limited number of spheres, where isolated spheres appear due to the algorithm looking for local maxima as opposed to global ones.
* **show_progress**: Print progress to the console.
* **initial_sphere_table**: Prior solver state. Can be in voxel units (if passed to `multisphere_from_voxels`) or physical units (if passed to `multisphere_from_mesh`). Must be a matrix of shape `(N, 4)` where each row is `(center_x, center_y, center_z, radius)`. If provided, the solver will start with this configuration instead of an empty table, which can significantly speed up convergence if you have a good initial guess.