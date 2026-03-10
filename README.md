<p align="center">
  <img src="https://raw.githubusercontent.com/FelixBuchele/multisphere/main/logo/multisphere_banner_ext.png"
       alt="multisphere logo"
       width="85%">
</p>

<p align="center">
  <a href="https://opensource.org/licenses/GPL-3.0">
    <img src="https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat-square"
         alt="License: GPLv3">
  </a>
  <a href="#citation">
    <img src="https://img.shields.io/badge/DOI-pending-lightgrey.svg?style=flat-square"
         alt="DOI: pending">
  </a>
</p>

---

# multisphere-cpp

`multisphere-cpp` creates overlapping-sphere representations of arbitrary 3D geometries based on voxelized Euclidean distance transforms (EDT) and feature-enhanced distance transform (FEDT) fields. The algorithm is designed for Discrete Element Method (DEM) simulations, where accurate yet computationally efficient particle shape representations are essential.

**This repository is a fork of the [Python implementation by Felix Buchele](https://github.com/FelixBuchele/multisphere).**

## Features

- Multisphere reconstruction from:
  - Triangle meshes (STL)
  - Binary voxel grids
- Exact EDT-driven sphere placement
- Iterative residual correction using FEDT
- Multiple termination criteria:
  - Shape precision
  - Maximum number of spheres
  - Minimum allowed sphere radius
- Robust Generalized Winding Number voxelization (via libigl)
- OpenMP parallelization for voxelization and field processing
- Export formats:
  - CSV (sphere centers & radii)
  - VTK (visualization)
  - STL (mesh export)

## Scientific Background

The multisphere algorithm is based on:
- Voxelization of the target geometry
- Exact Euclidean Distance Transform (EDT)
- Peak refinement of EDT maxima
- Iterative residual correction using the Feature-Enhanced Distance Transform (FEDT)
- Termination by shape accuracy, minimum radius, or maximum sphere count

The use of FEDT preserves the medial axis of the geometry and avoids the major drawbacks of greedy sphere removal methods, such as spurious small spheres, symmetry violations, and excessive runtime.

---

## C++ Implementation

The C++ implementation is designed for high-performance integration. It is a **header-only** library (core logic) with dependencies provided in `include/`.

### Dependencies

* **System**: CMake (≥3.10), C++17 compiler, OpenMP (optional but recommended).
* **Bundled (in `include/thirdparty/`)**: `libigl` (math/geometry), `cnpy` (numpy IO), `edt` (distance transform).
* **Not bundled:**  
  - [Eigen](https://eigen.tuxfamily.org/) (required, must be installed separately)
  - [zlib](https://zlib.net/) (required for `cnpy`, must be available on your system)

### Building the C++ Examples

Example scripts (`main.cpp`, `main_mesh.cpp`, etc.) are located in `src/` and can be built as follows:

```bash
mkdir build
cd build
cmake ../src
make -j4
```

> **Note:** Make sure `Eigen` and `zlib` are installed and discoverable by CMake before building.

### Using as a Header-Only Library



You can use `multisphere-cpp` as a header-only library in your own project:
- Add the `include/` directory to your compiler's include path.
- **Single include:** `#include "multisphere-interface.h"` gives you access to the entire public API.
- **All public API is in the `MSS` namespace.** You must either prefix all types and functions with `MSS::`, or add `using namespace MSS;` in your `.cpp` files.
- **Default argument values** for API functions are shown as comments in the interface header for clarity.
- **Note:** The `Eigen` library is required but **not provided** in the `include/` directory. You must have Eigen installed and available in your include path.
- **Note:** The `cnpy` library requires `zlib` to be available on your system.
- It is recommended to use the provided CMake configuration, or ensure your own CMake setup finds and links all required dependencies (`Eigen`, `zlib`, etc.) when including `multisphere-cpp` headers.
- No need to build the example executables unless you want to run the demos.


### Basic C++ Usage


The core API is provided via the umbrella header `multisphere-interface.h` and is in the `MSS` namespace.

```cpp
#include "multisphere-interface.h"

using namespace MSS; // Or use MSS:: prefix for all types/functions

int main() {
  // 1. Load Mesh
  FastMesh mesh = load_mesh_fast("example_mesh.stl");

  // 2. Run Reconstruction
  SpherePack sp = multisphere_from_mesh(
    mesh,
    150,    // div (resolution)
    2,      // padding
    8,      // min_radius_vox
    0.99,   // precision_target
    4,      // min_center_distance_vox
    100,    // max_spheres
    true,   // show_progress
    false   // confine_mesh
    // existing sphere pack can be passed as an additional argument
  );

  // 3. Export
  export_to_csv(sp, "results.csv");
  export_to_vtk(sp, "results.vtk");
  save_mesh_to_stl(grid_to_mesh(sp), "results.stl");
  return 0;
}
```

---

## Input/Output

- **Mesh loading:** STL files via `load_mesh_fast`
- **Voxel grid loading:** NumPy `.npy` files via `load_voxels_from_npy`
- **Export:** CSV, VTK, STL (no runtime visualization; use external tools for viewing)

### Visualization

You can visualize the output files using external tools such as [ParaView](https://www.paraview.org/) (for VTK), [MeshLab](https://www.meshlab.net/) (for STL), or any spreadsheet software (for CSV).  
No runtime visualization is included in this library.

---

## License

This project is licensed under the GNU General Public License v3.0.

See the LICENSE file for full details.

`multisphere-cpp` depends on third-party libraries with compatible licenses:

| Package      | License    | Usage         |
| ------------ | ---------- | ------------ |
| **Eigen**    | MPL2       | C++ Math     |
| **libigl**   | MPL2       | C++ Voxelization |
| **cnpy**     | MIT        | C++ IO       |
| **edt**      | MIT        | C++ Distance Transform |

## Author

**Arash Moradian**  
Friedrich-Alexander-Universität Erlangen–Nürnberg (FAU)  
Institute for Multiscale Simulation (MSS)  
moradian.arash@gmail.com

Contributors: Felix Buchele, Patric Müller, Thorsten Pöschel

## Citation

Felix Buchele, Patric Müller, Thorsten Pöschel,  
*Multi-Sphere-Shape generator for DEM simulations* manuscript in preparation

---

## Contact & Support

For questions, bug reports, or contributions, please [open an issue](https://github.com/aqa-arash/multisphere-cpp/issues) 