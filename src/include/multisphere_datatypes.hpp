
#ifndef MULTISPHERE_DATATYPES_HPP
#define MULTISPHERE_DATATYPES_HPP


/**
 * @file multisphere_datatypes.hpp
 * @brief Core data structures for multisphere-cpp library (namespace MSS).
 *
 * Defines FastMesh, SpherePack, VoxelGrid, and utility types for mesh and voxel operations.
 * Provides distance transform and sphere kernel generation.
 *
 * @author Arash Moradian
 * @date 2026-03-09
 */

#include <iostream>
#include <vector>
#include <thread>
#include <filesystem>
#ifdef HAVE_OPENMP
    #include <omp.h>
#endif
#include <memory>
#include <Eigen/Dense>
#include "thirdparty/edt.hpp" ///< High-performance C++ EDT library


namespace MSS {

/**
 * @brief Mesh structure for fast voxelization.
 */
struct FastMesh {
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

    // Global Physical Properties of the Multisphere Union
    float precision = 0.0f;           ///< Final precision achieved
    float volume = 0.0f;
    Eigen::Vector3f center_of_mass = Eigen::Vector3f::Zero();
    Eigen::Matrix3f inertia_tensor = Eigen::Matrix3f::Zero();
    Eigen::Matrix3f principal_axes = Eigen::Matrix3f::Identity();
    Eigen::Vector3f principal_moments = Eigen::Vector3f::Zero();

    SpherePack(Eigen::MatrixX3f c = Eigen::MatrixX3f(0, 3), Eigen::VectorXf r = Eigen::VectorXf(0)) 
        : centers(std::move(c)), radii(std::move(r)) {
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

    inline typename std::vector<T>::reference operator()(size_t x, size_t y, size_t z) {
        return data[x * (shape[1] * shape[2]) + y * shape[2] + z];
    }

    inline typename std::vector<T>::const_reference operator()(size_t x, size_t y, size_t z) const {
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

        std::unique_ptr<float[]> dists = nullptr;

        if (binary_mode) {
            // Binary mode: Flatten all labels > 0 into a uniform mask of 1s
            std::vector<uint8_t> temp_mask(this->data.size());
            #pragma omp parallel for
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
        } else {
            // skip for now, as we currently only use binary mode. Future: implement multi-label EDT if needed.
            throw std::runtime_error("Non-binary EDT mode not implemented yet.");
        }
        
        #pragma omp parallel for
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

/**
 * @brief Vertex key for hashing.
 */
struct VertexKey {
    float x, y, z;
    bool operator==(const VertexKey& other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

/**
 * @brief Hash function for VertexKey.
 */
struct VertexKeyHash {
    std::size_t operator()(const VertexKey& k) const {
        return std::hash<float>()(k.x) ^ 
              (std::hash<float>()(k.y) << 1) ^ 
              (std::hash<float>()(k.z) << 2);
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

#endif // MULTISPHERE_DATATYPES_HPP