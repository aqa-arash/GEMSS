#ifndef MULTISPHERE_DATATYPES_HPP
#define MULTISPHERE_DATATYPES_HPP

#include <iostream>
#include <vector>
#include <thread>
#ifdef HAVE_OPENMP
    #include <omp.h>
#endif
#include <memory>
#include <Eigen/Dense>
#include "include/edt.hpp" // High-performance C++ EDT library


struct FastMesh {
    // We use float for vertices as STL precision is usually 32-bit
    Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> vertices;
    Eigen::Matrix<int,   Eigen::Dynamic, 3, Eigen::RowMajor> triangles;
    // Check if the mesh is valid for voxelization
    bool is_empty() const { return vertices.size() == 0; }
};


// --- SpherePack ---
struct SpherePack {
    Eigen::MatrixX3f centers;
    Eigen::VectorXf radii;

    SpherePack(Eigen::MatrixX3f c = Eigen::MatrixX3f(0, 3), Eigen::VectorXf r = Eigen::VectorXf(0)) 
        : centers(std::move(c)), radii(std::move(r)) {
        if (centers.rows() != radii.size()) {
            throw std::invalid_argument("Centers and radii length mismatch.");
        }
    }

    // Modern C++ version of @property
    size_t num_spheres() const { return radii.size(); }
    float min_radius() const { return radii.size() > 0 ? radii.minCoeff() : 0.0; }
    float max_radius() const { return radii.size() > 0 ? radii.maxCoeff() : 0.0; }
};

// --- VoxelGrid ---
template<typename T>
class VoxelGrid {
public:
    std::vector<T> data;
    std::array<size_t, 3> shape; // {nx, ny, nz}
    float voxel_size;
    Eigen::Vector3f origin;

    VoxelGrid(size_t nx, size_t ny, size_t nz, float v_size = 1.0f, Eigen::Vector3f orig = Eigen::Vector3f::Zero())
        : shape({nx, ny, nz}), voxel_size(v_size), origin(orig) {
        data.resize(nx * ny * nz, static_cast<T>(0));
    }

    inline typename std::vector<T>::reference operator()(size_t x, size_t y, size_t z) {
        return data[x * (shape[1] * shape[2]) + y * shape[2] + z];
    }

    inline typename std::vector<T>::const_reference operator()(size_t x, size_t y, size_t z) const {
        return data[x * (shape[1] * shape[2]) + y * shape[2] + z];
    }

    size_t nx() const { return shape[0]; }
    size_t ny() const { return shape[1]; }
    size_t nz() const { return shape[2]; }

    // Distance Transform using the 'edt' library
    // This returns a VoxelGrid of floats containing physical distances
   // Distance Transform using the 'edt' library
    VoxelGrid<float> distance_transform() const {
        VoxelGrid<float> result(nx(), ny(), nz(), voxel_size, origin);

        const uint8_t* data_ptr = nullptr;
        std::vector<uint8_t> temp_mask; // Keep temporary buffer alive if needed

        // OPTIMIZATION: Zero-copy if we are already a uint8_t grid
        // We use 'if constexpr' so this compiles cleanly even for float grids
        if constexpr (std::is_same<T, uint8_t>::value) {
            // Direct pointer access - NO COPY, NO ALLOCATION
            data_ptr = reinterpret_cast<const uint8_t*>(this->data.data());
        } 
        else {
            // Fallback for float/bool grids: Convert to uint8_t buffer
            temp_mask.resize(this->data.size());
            
            #pragma omp parallel for
            for(size_t i = 0; i < this->data.size(); ++i) {
                // Determine threshold: 0 vs >0
                temp_mask[i] = (this->data[i] > static_cast<T>(0)) ? 1 : 0;
            }
            data_ptr = temp_mask.data();
        }

        // Setup threads
        int num_threads = 1;
        #ifdef HAVE_OPENMP
            num_threads = omp_get_max_threads();
        #else
            num_threads = std::thread::hardware_concurrency();
            if (num_threads == 0) num_threads = 1;
        #endif

        // Call EDT library
        // Note: shape indices [2], [1], [0] (z, y, x) match the library's stride expectation
        float* dists = edt::edt<uint8_t>(
            const_cast<uint8_t*>(data_ptr), 
            static_cast<int>(shape[2]), 
            static_cast<int>(shape[1]), 
            static_cast<int>(shape[0]), 
            1.0f, 1.0f, 1.0f, 
            true, num_threads, nullptr
        );
        
        // Copy back to result (Float grid)
        // We cannot avoid this copy as 'edt' mallocs its own float buffer
        #pragma omp parallel for
        for(size_t i = 0; i < this->data.size(); ++i) {
            result.data[i] = dists[i];
        }

        free(dists); // Library uses malloc, so we must use free()
        return result;
    }

    // Factory: Sphere Kernel (Vectorized-style loop)
   // --- Static Factory: Sphere Kernel ---
    // Templated on <K> to allow creating kernels of bool, float, etc.
    // Usage: VoxelGrid<float>::sphere_kernel<float>(5, 1.0);
    void sphere_kernel(float cx, float cy, float cz, float radius, T fill_value = static_cast<T>(1)) {
        if (radius <= 0) return;

        float r_sq = radius * radius;

        // 1. Calculate Bounding Box (Clamped to Grid Dimensions)
        // We use member variables directly (nx_, ny_, nz_ or nx(), ny(), nz())
        // Assuming accessors nx(), ny(), nz() exist based on your file structure
        int n_x = static_cast<int>(this->nx());
        int n_y = static_cast<int>(this->ny());
        int n_z = static_cast<int>(this->nz());

        int min_x = std::max(0, static_cast<int>(std::floor(cx - radius)));
        int max_x = std::min(n_x - 1, static_cast<int>(std::ceil(cx + radius)));
        
        int min_y = std::max(0, static_cast<int>(std::floor(cy - radius)));
        int max_y = std::min(n_y - 1, static_cast<int>(std::ceil(cy + radius)));
        
        int min_z = std::max(0, static_cast<int>(std::floor(cz - radius)));
        int max_z = std::min(n_z - 1, static_cast<int>(std::ceil(cz + radius)));

        // 2. Optimized Rasterization Loop
        for (int x = min_x; x <= max_x; ++x) {
            float dx = x - cx;
            float dx2 = dx * dx;
            if (dx2 > r_sq) continue;

            for (int y = min_y; y <= max_y; ++y) {
                float dy = y - cy;
                float dy2 = dy * dy;
                float dxy2 = dx2 + dy2;
                if (dxy2 > r_sq) continue;

                for (int z = min_z; z <= max_z; ++z) {
                    float dz = z - cz;
                    float dist_sq = dxy2 + (dz * dz);

                    if (dist_sq <= r_sq) {
                        // Direct member access. 
                        // Assuming operator() returns a reference to the data.
                        (*this)(x, y, z) = fill_value;
                    }
                }
            }
        }
    }
};


// Custom Hash for Vertex Key
struct VertexKey {
    float x, y, z;
    
    bool operator==(const VertexKey& other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

// Hasher for unordered_map
struct VertexKeyHash {
    std::size_t operator()(const VertexKey& k) const {
        // Simple XOR hash
        return std::hash<float>()(k.x) ^ 
              (std::hash<float>()(k.y) << 1) ^ 
              (std::hash<float>()(k.z) << 2);
    }
};

#endif