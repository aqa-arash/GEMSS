<p align="center">
  <img src="https://raw.githubusercontent.com/aqa-arash/GEMSS/463b63bf356f5c7a27b46df3f7771c4d2510eb00/logo/multisphere_banner_ext.png"
       alt="multisphere logo"
       width="85%">
</p>

<p align="center">
  <a href="https://opensource.org/licenses/MIT">
    <img src="https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square"
         alt="License: MIT">
  </a>
  <a href="https://doi.org/10.1016/j.powtec.2026.123162">
    <img src="https://img.shields.io/badge/DOI-10.1016%2Fj.powtec.2026.123162-blue.svg?style=flat-square" alt="DOI">
  </a>
</p>
---

# GEMSS - GEnerator of Multi-Sphere Shapes

`GEMSS` is a header only C++17 library that creates overlapping-sphere representations of arbitrary 3D geometries based on voxelized Euclidean distance transforms (EDT) and feature-enhanced distance transform (FEDT) fields. It equips the output with physical properties—including Center of Mass (CoM), Inertia Tensors, and Principal Moments—necessary for Discrete Element Method (DEM) and Molecular Dynamics (MD) simulations.


## Features

- Multisphere reconstruction according to MSS algorithm [1] from:
  - Triangle meshes (binary STL) 
  - Binary voxel grids
- Physical Property Computation: Automated calculation of volume, Center of Mass (CoM), inertia tensors, and principal axes directly from the multisphere union or the target geometry, along with principal extents (L ≥ I ≥ S) calculated from the reconstructed sphere geometry.
- Topology Filtering: Built-in pruning of isolated sphere networks to guarantee continuous rigid-body representations.
- Multiple termination criteria:
  - Shape precision
  - Maximum number of spheres
  - Minimum allowed sphere radius
- Robust Generalized Winding Number voxelization (via libigl)
- OpenMP parallelization for voxelization and field processing
- Export formats:
  - JSON (SpherePack)
  - CSV (sphere centers & radii)
  - VTK (visualization)
  - STL (mesh export)


## Scientific Background

The MSS algorithm is based on:
- Voxelization of the target geometry
- Euclidean Distance Transform (EDT)
- Iterative residual correction using the Feature-Enhanced Distance Transform (FEDT)
- Termination by shape accuracy, minimum radius, or maximum sphere count

The use of FEDT preserves the medial axis of the geometry and avoids the major drawbacks of greedy sphere removal methods, such as spurious small spheres and symmetry violations [2].

---

## C++ Implementation

The C++ implementation is designed for high-performance integration. It is a **header-only** library. All necessary third-party mathematics and geometry processing headers are bundled in the `thirdparty/` directory, with the sole exception of Eigen.

### ⚠️ Critical Performance Requirement

Because `GEMSS` relies heavily on massive voxel-space iterations and 3D distance transforms, **you should compile your code with aggressive optimizations enabled (`-O3` on GCC/Clang, or `/Ox` on MSVC)**. 

Without optimizations, the compiler will not inline the geometry math or vectorize the spatial hashing loops. Processing 1,000,000 query points takes approximately **10 seconds** with `-O3`, but will take **over 10 minutes** without it. Do not evaluate the performance of this library in an unoptimized state.

### Dependencies

* **System**: CMake (≥3.15), C++17 compiler, OpenMP (optional but highly recommended).
* **Bundled (in `GEMSS/thirdparty/`)**: `libigl` (math/geometry), `edt` (distance transform).
* **Not bundled:** [Eigen](https://eigen.tuxfamily.org/) (required, must be installed separately).

---

## Step-by-Step Integration Guide

Because this is a header-only library, you do not need to pre-compile anything to use it in your own software. Below is the exact process to integrate and compile your own code using `GEMSS`.

### Option A: Using CMake (Recommended)

**Step 1:** Set up your project directory. Copy the `GEMSS` folder into your project's `include` directory (or any directory you use for external headers).
```text
my_project/
├── CMakeLists.txt
├── main.cpp
└── include/
    └── GEMSS/  <-- Copy GEMSS here
```

**Step 2:** Update your `CMakeLists.txt` to add the GEMSS directory to your include paths and link Eigen.
```cmake
cmake_minimum_required(VERSION 3.15)
project(MyMultisphereApp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required external dependency
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Add GEMSS header-only library
add_subdirectory(include/GEMSS)

# Create your executable(s)
add_executable(my_app main.cpp)

# Link GEMSS header-only library
# GEMSS handles all include paths internally
# Link GEMSS (includes Eigen and configured OpenMP dependencies)
target_link_libraries(my_app PRIVATE GEMSS::GEMSS)

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

If you do not want to use CMake, you can compile your code directly via the terminal. You must explicitly provide the include paths for GEMSS, its bundled `thirdparty` files, and `Eigen3`.

Assuming GEMSS is in your include directory:

```bash
g++ -std=c++17 -O3 -march=native -fopenmp main.cpp \
    -I include/GEMSS/ \
    -I include/GEMSS/thirdparty/ \
    -I include/GEMSS/thirdparty/igl \
    -I /usr/include/eigen3/ \
    -o my_app
```
*(Note: Replace `/usr/include/eigen3/` with the actual path to your Eigen installation if it differs).*

#### Enabling Internal Debug Prints

If you need to track internal library execution (such as function calls and intermediate values), enable the `GEMSS_DEBUG` flag. 

If using CMake, pass the option during configuration:
```bash
cmake ../src -DGEMSS_DEBUG=ON
```

If compiling directly with G++, pass the definition flag:
```bash
g++ -std=c++17 -O3 -march=native -fopenmp -DGEMSS_DEBUG main.cpp ...
```

---

## Basic C++ Usage

The core API is provided via the umbrella header `GEMSS-interface.h` and is contained entirely within the `GEMSS` namespace.

```cpp
#include "GEMSS-interface.h"
#include <iostream>

using namespace GEMSS;

int main() {
  // 1. Load Mesh
  STLMesh mesh = load_mesh("example_mesh.stl");

  // 2. Set up configuration
  MultisphereConfig config;
  config.div = 150;                   // Voxel grid resolution
  config.padding = 2;                 // Grid padding
  config.min_radius_vox = 8;          // Minimum sphere radius in voxels
  config.precision_target = 0.99;     // Target precision
  config.min_center_distance_rel = 1.4f; // Relative sphere center spacing
  config.max_spheres = 100;           // Maximum number of spheres
  config.show_progress = true;        // Show progress output
  config.confine_mesh = false;        // Do not confine spheres to mesh boundary
  config.prune_isolated_spheres = true; // Remove disconnected spheres
  config.compute_physics = 1;         // Compute physics (1 = based on reconstruction)

  // 3. Run Reconstruction
  SpherePack sp = multisphere_from_mesh(mesh, config);

  // 4. Access Computed Physics Properties
  std::cout << "Mass: " << sp.mass << "\n";
  std::cout << "Center of Mass:\n" << sp.center_of_mass << "\n";
  std::cout << "Inertia Tensor:\n" << sp.inertia_tensor << "\n";
  std::cout << "Principal Extents (L, I, S): "
          << sp.principal_extents.transpose() << "\n";

  // 5. Export
  export_to_csv(sp, "results.csv");
  export_to_vtk(sp, "results.vtk");

  return 0;
}
```

---

## Input/Output

- **Mesh loading:** STL (binary) files via `load_mesh`
- **Export:** JSON, CSV, VTK, STL (no runtime visualization; use external tools for viewing)

### Visualization

You can visualize the output files using external tools such as [ParaView](https://www.paraview.org/) (for VTK), [MeshLab](https://www.meshlab.net/) (for STL), or any spreadsheet software (for CSV).  
No runtime visualization is included in this library.


### Documentation

- [Configuration](./README_CONFIG.md) — Configuration parameters and defaults.
- [Data Types](./README_DATATYPES.md) — GEMSS data structures, including `principal_extents`.
- [Examples](./README_EXAMPLES.md) — C++ and Python examples.

---

## License

The core of this project is licensed under the MIT License. However, it includes third-party components located in the third_party/ directory that are subject to the MPLv2 license.

See the `LICENSE` file for full details.

`multisphere-cpp` depends on third-party libraries with compatible licenses:

| Package | License | Usage |
| :--- | :--- | :--- |
| **Eigen** | MPL2 | C++ Math |
| **libigl** | MPL2 | C++ Voxelization |


## Author

**Arash Moradian** Friedrich-Alexander-Universität Erlangen–Nürnberg (FAU)  
Institute for Multiscale Simulation (MSS)  
moradian.arash@gmail.com  

Contributors: Felix Buchele, Patric Müller, Thorsten Pöschel


## Citation

Algorithm :
Felix Buchele, Thorsten Pöschel, Patric Müller. *Multi-sphere shape generator for DEM simulations of complex-shaped particles*. Powder Technology, 480, 122635 (2026). https://doi.org/10.1016/j.powtec.2026.122635

Proof of Medial Axis:
Arash Moradian, Felix Buchele, Thorsten Poeschel. Why the Multi-Sphere Shape Generator Works: Medial-Axis Placement of Spheres. arXiv preprint arXiv:2608.15457 (2026). https://arxiv.org/abs/2608.15457

Implementation : 
A. Moradian, F. Buchele, P. Müller, T. Pöschel. *GEMSS: A C++ library for multi-sphere modeling in DEM simulations*. Powder Technology, 123162 (2026). https://doi.org/10.1016/j.powtec.2026.123162

---

## Contact & Support

For questions, bug reports, or contributions, please [open an issue](https://github.com/aqa-arash/GEMSS/issues)
