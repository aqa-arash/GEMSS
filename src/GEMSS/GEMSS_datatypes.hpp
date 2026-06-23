/**
 * @file GEMSS_datatypes.hpp
 * @brief Core data structures for multisphere-cpp library (namespace GEMSS).
 *
 * Defines STLMesh, SpherePack, VoxelGrid, and utility types for mesh and voxel operations.
 * Provides distance transform and sphere kernel generation.
 *
 * @author Arash Moradian
 * @date 2026-03-09
 */

#ifndef GEMSS_DATATYPES_HPP
#define GEMSS_DATATYPES_HPP

#include <iostream>
#include <vector>
#include <thread>
#include <filesystem>
#ifdef HAVE_OPENMP
    #include <omp.h>
#endif
#include <memory>
#include <Eigen/Core>
#include "thirdparty/edt.hpp" ///< High-performance C++ EDT library


namespace GEMSS {

/**
 * @brief Mesh structures.
 */
struct STLMesh {
    Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> vertices; ///< Mesh vertices (float, 3D)
    Eigen::Matrix<int,   Eigen::Dynamic, 3, Eigen::RowMajor> triangles; ///< Mesh triangles (int, 3D)
    /**
     * @brief Check if the mesh is empty.
     * @return True if mesh has no vertices.
     */
    bool is_empty() const { return vertices.size() == 0; }
};

/**
 * @brief Sphere pack structure for reconstruction results.
 */
struct SpherePack {
    Eigen::MatrixX3f centers; ///< Sphere centers (Nx3)
    Eigen::VectorXf radii;    ///< Sphere radii (N)
    bool is_2d = false; 

    bool is_2d = false;   // planar pack: physics is area/Izz, split collapses out-of-plane axis

    // Global Physical Properties of the Multisphere Union
    float precision = 0.0f;           ///< Final precision achieved
    float density = 1.0f;             ///< Material density for mass properties
    float mass = 0.0f;                ///< Mass of the union of spheres
    float bounding_radius= 0.0f;  // Bounding sphere radius of the union
    Eigen::Vector3f center_of_mass = Eigen::Vector3f::Zero();
    Eigen::Matrix3f inertia_tensor = Eigen::Matrix3f::Zero();
    Eigen::Matrix3f principal_axes = Eigen::Matrix3f::Identity();
    Eigen::Vector3f principal_moments = Eigen::Vector3f::Zero();

    SpherePack(Eigen::MatrixX3f c = Eigen::MatrixX3f(0, 3), Eigen::VectorXf r = Eigen::VectorXf(0)) 
        : centers (c), radii(r) {
        if (centers.rows() != radii.size()) {
            throw std::invalid_argument("Centers and radii length mismatch.");
        }
    }

    size_t num_spheres() const { return radii.size(); }
    float min_radius() const { return radii.size() > 0 ? radii.minCoeff() : 0.0f; }
    float max_radius() const { return radii.size() > 0 ? radii.maxCoeff() : 0.0f; }
};

/**
 * @brief Voxel grid template for 3D data.
 * @tparam T Data type (bool, float, uint8_t, etc.)
 */
template<typename T>
class VoxelGrid {
public:
    std::vector<T> data;              ///< Voxel data buffer
    std::array<size_t, 3> shape;      ///< Grid shape {nx, ny, nz}
    float voxel_size;                 ///< Physical voxel size
    Eigen::Vector3f origin;           ///< Grid origin

    VoxelGrid(size_t nx, size_t ny, size_t nz, float v_size = 1.0f, Eigen::Vector3f orig = Eigen::Vector3f::Zero())
        : shape({nx, ny, nz}), voxel_size(v_size), origin(orig) {
        data.resize(nx * ny * nz, static_cast<T>(0));
    }

    /**
     * @brief Named constructor for a 2D grid (single z-slice, nz == 1).
     * A plain 2-arg ctor would be ambiguous with the 3D ctor, so use this.
     */
    static VoxelGrid<T> make_2d(size_t nx, size_t ny,
                                float v_size = 1.0f,
                                Eigen::Vector3f orig = Eigen::Vector3f::Zero()) {
        return VoxelGrid<T>(nx, ny, 1, v_size, orig);
    }

    /// @brief True if this grid is a single z-slice (treated as 2D).
    bool is_2d() const { return shape[2] == 1; }

    inline typename std::vector<T>::reference operator()(size_t x, size_t y, size_t z = 0) {
        return data[x * (shape[1] * shape[2]) + y * shape[2] + z];
    }

    inline typename std::vector<T>::const_reference operator()(size_t x, size_t y, size_t z = 0) const {
        return data[x * (shape[1] * shape[2]) + y * shape[2] + z];
    }

    size_t nx() const { return shape[0]; }
    size_t ny() const { return shape[1]; }
    size_t nz() const { return shape[2]; }

    /**
     * @brief Compute distance transform using EDT library.
     * @return VoxelGrid<float> with physical distances.
     */
    VoxelGrid<float> distance_transform(bool binary_mode = true ) const {
        VoxelGrid<float> result(nx(), ny(), nz(), voxel_size, origin);

        int num_threads = 1;
        #ifdef HAVE_OPENMP
            num_threads = omp_get_max_threads();
        #else
            num_threads = std::thread::hardware_concurrency();
            if (num_threads == 0) num_threads = 1;
        #endif

        // --- 2D fast path: a single-slice grid must use the 2D EDT. ---
        // The 3D EDT with black_border on a 1-voxel-thick slab would clamp every
        // foreground voxel to distance 1 from the empty z-faces, destroying the field.
        if (shape[2] == 1) {
            std::vector<uint8_t> mask(this->data.size());
            for (size_t i = 0; i < this->data.size(); ++i)
                mask[i] = (this->data[i] > static_cast<T>(0)) ? 1 : 0;

            // Layout idx = x*ny + y  =>  y is contiguous. EDT wants x fastest,
            // so EDT sx = ny (fast axis), EDT sy = nx.
            std::unique_ptr<float[]> dists2d(edt::binary_edt<uint8_t>(
                mask.data(),
                static_cast<int>(shape[1]),   // sx (fast)
                static_cast<int>(shape[0]),   // sy
                1.0f, 1.0f,
                /*black_border=*/true, num_threads, nullptr));

            for (size_t i = 0; i < this->data.size(); ++i) result.data[i] = dists2d[i];
            return result;
        }
        
        std::unique_ptr<float[]> dists = nullptr;

        if (binary_mode) {
            if constexpr (std::is_same_v<T, uint8_t>) {
                // 1. Cross-platform OpenMP check for multi-label data
                int needs_flattening = 0;
                
                #pragma omp parallel for reduction(|:needs_flattening)
                for (size_t i = 0; i < this->data.size(); ++i) {
                    if (this->data[i] > 1) {
                        needs_flattening |= 1;
                    }
                }

                // 2. Branch: Flatten multi-label data OR use zero-copy fast path
                if (needs_flattening != 0) {
                    #ifdef MULTISPHERE_DEBUG
                    std::cerr << "[Warning] VoxelGrid uint8_t data contains multi-label values (> 1). Flattening to strict binary mask." << std::endl;
                    #endif
                    
                    std::vector<uint8_t> temp_mask(this->data.size());
                    for(size_t i = 0; i < this->data.size(); ++i) {
                        temp_mask[i] = (this->data[i] > 0) ? 1 : 0;
                    }
                    
                    dists.reset(edt::binary_edt<uint8_t>(
                        temp_mask.data(), 
                        static_cast<int>(shape[2]), 
                        static_cast<int>(shape[1]), 
                        static_cast<int>(shape[0]), 
                        1.0f, 1.0f, 1.0f, 
                        true, num_threads, nullptr
                    ));
                } else {
                    // Zero-copy fast path
                    dists.reset(edt::binary_edt<uint8_t>(
                        const_cast<uint8_t*>(this->data.data()), 
                        static_cast<int>(shape[2]), 
                        static_cast<int>(shape[1]), 
                        static_cast<int>(shape[0]), 
                        1.0f, 1.0f, 1.0f, 
                        true, num_threads, nullptr
                    ));
                }
            } else {
                // Safe allocation block for generic non-uint8_t types
                std::vector<uint8_t> temp_mask(this->data.size());
                for(size_t i = 0; i < this->data.size(); ++i) {
                    temp_mask[i] = (this->data[i] > static_cast<T>(0)) ? 1 : 0;
                }
                
                dists.reset(edt::binary_edt<uint8_t>(
                    temp_mask.data(), 
                    static_cast<int>(shape[2]), 
                    static_cast<int>(shape[1]), 
                    static_cast<int>(shape[0]), 
                    1.0f, 1.0f, 1.0f, 
                    true, num_threads, nullptr
                ));
            }
        } else {
            throw std::runtime_error("Non-binary EDT mode not implemented yet.");
        }
        
        
        for(size_t i = 0; i < this->data.size(); ++i) {
            result.data[i] = dists[i];
        }

        return result;
    }
    

    /**
     * @brief Fill grid with a sphere kernel.
     * @param cx Center x
     * @param cy Center y
     * @param cz Center z
     * @param radius Sphere radius
     * @param fill_value Value to fill inside sphere
     */
    void sphere_kernel(float cx, float cy, float cz, float radius, T fill_value = static_cast<T>(1)) {
        if (radius <= 0) return;

        float r_sq = radius * radius;
        int n_x = static_cast<int>(this->nx());
        int n_y = static_cast<int>(this->ny());
        int n_z = static_cast<int>(this->nz());

        int min_x = std::max(0, static_cast<int>(std::floor(cx - radius)));
        int max_x = std::min(n_x - 1, static_cast<int>(std::ceil(cx + radius)));
        int min_y = std::max(0, static_cast<int>(std::floor(cy - radius)));
        int max_y = std::min(n_y - 1, static_cast<int>(std::ceil(cy + radius)));
        int min_z = std::max(0, static_cast<int>(std::floor(cz - radius)));
        int max_z = std::min(n_z - 1, static_cast<int>(std::ceil(cz + radius)));

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
                        (*this)(x, y, z) = fill_value;
                    }
                }
            }
        }
    }
};


// Helper struct for the Sort-and-Sweep algorithm
struct RawVertex {
    float x, y, z;
    uint32_t original_id;

    // Strict weak ordering for std::sort
    bool operator<(const RawVertex& other) const {
        if (x != other.x) return x < other.x;
        if (y != other.y) return y < other.y;
        return z < other.z;
    }

    // Exact equality check
    bool operator==(const RawVertex& other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

} // namespace MSS

#endif // GEMSS_DATATYPES_HPP
