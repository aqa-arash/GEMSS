<p align="center">
  <img src="[https://raw.githubusercontent.com/FelixBuchele/multisphere/main/logo/multisphere_banner_ext.png](https://raw.githubusercontent.com/FelixBuchele/multisphere/main/logo/multisphere_banner_ext.png)"
       alt="multisphere logo"
       width="85%">
</p>

<p align="center">
  <a href="[https://opensource.org/licenses/GPL-3.0](https://opensource.org/licenses/GPL-3.0)">
    <img src="[https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat-square](https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat-square)"
         alt="License: GPLv3">
  </a>
  <a href="#citation">
    <img src="[https://img.shields.io/badge/DOI-pending-lightgrey.svg?style=flat-square](https://img.shields.io/badge/DOI-pending-lightgrey.svg?style=flat-square)"
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
- **Physical Property Computation:** Automated calculation of volume, Center of Mass (CoM), inertia tensors, and principal axes directly from the multisphere union or the target geometry.
- **Topology Filtering:** Built-in pruning of isolated sphere networks to guarantee continuous rigid-body representations.
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

The C++ implementation is designed for high-performance integration. It is a **pure header-only** library. All necessary third-party mathematics and geometry processing headers are bundled in the `include/` directory, with the sole exception of Eigen.

### ⚠️ Critical Performance Requirement

Because `multisphere-cpp` relies heavily on massive voxel-space iterations and 3D distance transforms, **you must compile your code with aggressive optimizations enabled (`-O3` on GCC/Clang, or `/Ox` on MSVC)**. 

Without optimizations, the compiler will not inline the geometry math or vectorize the spatial hashing loops. Processing 1,000,000 query points takes approximately **10 seconds** with `-O3`, but will take **over 10 minutes** without it. Do not evaluate the performance of this library in an unoptimized state.

### Dependencies

* **System**: CMake (≥3.15), C++17 compiler, OpenMP (optional but highly recommended).
* **Bundled (in `include/thirdparty/`)**: `libigl` (math/geometry), `edt` (distance transform).
* **Not bundled:** [Eigen](https://eigen.tuxfamily.org/) (required, must be installed separately).

---

## Step-by-Step Integration Guide

Because this is a header-only library, you do not need to pre-compile anything to use it in your own software. Below is the exact process to integrate and compile your own code using `multisphere-cpp`.

### Option A: Using CMake (Recommended)

**Step 1:** Set up your project directory. Clone this repository into a `thirdparty` folder inside your project.
```text
my_project/
├── CMakeLists.txt
├── main.cpp
└── thirdparty/
    └── multisphere-cpp/  <-- Cloned repository here
```

**Step 2:** Write your `CMakeLists.txt`. This file explicitly requests `Eigen3` and links the library.
```cmake
cmake_minimum_required(VERSION 3.15)
project(MyMultisphereApp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required external dependency
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Add the multisphere-cpp directory (disabling its example builds)
set(MULTISPHERE_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory(thirdparty/multisphere-cpp)

# Create your executable
add_executable(my_app main.cpp)

# Link the header-only library
target_link_libraries(my_app PRIVATE multisphere_lib)

# FORCE OPTIMIZATIONS for your executable (Critical!)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(my_app PRIVATE -O3 -march=native)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(my_app PRIVATE /Ox /Oi /Ot)
endif()
```

**Step 3:** Build your code from the terminal.
```bash
cd my_project
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j4
./my_app
```

### Option B: Compiling Directly with G++ (Command Line)

If you do not want to use CMake, you can compile your code directly via the terminal. You must explicitly provide the include paths for `multisphere-cpp`, its bundled `thirdparty` files, and `Eigen3`.

Assuming `multisphere-cpp` is in the same directory as your `main.cpp`, run the following exact command:

```bash
g++ -std=c++17 -O3 -march=native -fopenmp main.cpp \
    -I multisphere-cpp/ \
    -I multisphere-cpp/include/ \
    -I multisphere-cpp/include/thirdparty/ \
    -I multisphere-cpp/include/thirdparty/igl \
    -I /usr/include/eigen3/ \
    -o my_app
```
*(Note: Replace `/usr/include/eigen3/` with the actual path to your Eigen installation if it differs).*

#### Enabling Internal Debug Prints

If you need to track internal library execution (such as function calls and intermediate values), enable the `MULTISPHERE_DEBUG` flag. 

If using CMake, pass the option during configuration:
```bash
cmake .. -DCMAKE_BUILD_TYPE=Release -DMULTISPHERE_DEBUG=ON
```

If compiling directly with G++, pass the definition flag:
```bash
g++ -std=c++17 -O3 -march=native -fopenmp -DMULTISPHERE_DEBUG main.cpp ...
```

---

## Basic C++ Usage

The core API is provided via the umbrella header `multisphere-interface.h` and is contained entirely within the `MSS` namespace.

```cpp
#include "multisphere-interface.h"
#include <iostream>

using namespace MSS;

int main() {
  // 1. Load Mesh
  FastMesh mesh = load_mesh_fast("example_mesh.stl");

  // 2. Set up configuration
  MultisphereConfig config;
  config.div = 150;                   // Voxel grid resolution
  config.padding = 2;                 // Grid padding
  config.min_radius_vox = 8;          // Minimum sphere radius in voxels
  config.precision_target = 0.99;     // Target precision
  config.min_center_distance_vox = 4; // Minimum center distance in voxels
  config.max_spheres = 100;           // Maximum number of spheres
  config.show_progress = true;        // Show progress output
  config.confine_mesh = false;        // Do not confine spheres to mesh boundary
  config.prune_isolated_spheres = true; // Remove disconnected spheres
  config.compute_physics = 1;         // Compute physics (1 = based on reconstruction)

  // 3. Run Reconstruction
  SpherePack sp = multisphere_from_mesh(mesh, config);

  // 4. Access Computed Physics Properties
  std::cout << "Reconstruction Volume: " << sp.volume << "\n";
  std::cout << "Center of Mass:\n" << sp.center_of_mass << "\n";
  std::cout << "Inertia Tensor:\n" << sp.inertia_tensor << "\n";

  // 5. Export
  export_to_csv(sp, "results.csv");
  export_to_vtk(sp, "results.vtk");

  return 0;
}
```

---

## Input/Output

- **Mesh loading:** STL files via `load_mesh_fast`
- **Export:** CSV, VTK, STL (no runtime visualization; use external tools for viewing)

### Visualization

You can visualize the output files using external tools such as [ParaView](https://www.paraview.org/) (for VTK), [MeshLab](https://www.meshlab.net/) (for STL), or any spreadsheet software (for CSV).  
No runtime visualization is included in this library.

---

## License

This project is licensed under the GNU General Public License v3.0.

See the `LICENSE` file for full details.

`multisphere-cpp` depends on third-party libraries with compatible licenses:

| Package | License | Usage |
| :--- | :--- | :--- |
| **Eigen** | MPL2 | C++ Math |
| **libigl** | MPL2 | C++ Voxelization |
| **edt** | MIT | C++ Distance Transform |

## Author

**Arash Moradian** Friedrich-Alexander-Universität Erlangen–Nürnberg (FAU)  
Institute for Multiscale Simulation (MSS)  
moradian.arash@gmail.com  

Contributors: Felix Buchele, Patric Müller, Thorsten Pöschel

## Citation

Felix Buchele, Patric Müller, Thorsten Pöschel,  
*Multi-Sphere-Shape generator for DEM simulations* manuscript in preparation

---

## Contact & Support

For questions, bug reports, or contributions, please [open an issue](https://github.com/aqa-arash/multisphere-cpp/issues)