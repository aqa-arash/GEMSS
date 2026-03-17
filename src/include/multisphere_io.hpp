 #ifndef MULTISPHERE_IO_HPP
 #define MULTISPHERE_IO_HPP

/**
 * @file multisphere_io.hpp
 * @brief Input/output utilities for multisphere-cpp library.
 *
 * Provides mesh loading, sphere pack export, voxel grid export, and STL/NPY file utilities.
 *
 * @author Arash Moradian
 * @date 2026-03-09
 */

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <Eigen/Dense>
#include "thirdparty/cnpy.h"

#include "multisphere_datatypes.hpp"

namespace MSS {

/**
 * @brief Loads a mesh from a binary STL file.
 * @param path Path to STL file.
 * @return FastMesh structure.
 */
inline FastMesh load_mesh_fast(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) throw std::runtime_error("File not found: " + path);

    char header_check[6] = {0};
    file.read(header_check, 5);
    if (std::string(header_check) == "solid") {
        throw std::runtime_error("ASCII STL not supported. Please provide a Binary STL.");    
    }

    file.seekg(80);
    uint32_t num_triangles;
    file.read(reinterpret_cast<char*>(&num_triangles), 4);

    // O(1) file size validation
    std::error_code ec;
    uintmax_t physical_size = std::filesystem::file_size(path, ec);
    if (!ec && physical_size < 84 + (static_cast<uintmax_t>(num_triangles) * 50)) {
        throw std::runtime_error("Malformed Binary STL: File too small for claimed triangle count.");
    }

    if (num_triangles == 0) return FastMesh();

    const uint32_t total_raw_vertices = num_triangles * 3;
    std::vector<RawVertex> raw_verts;
    
    try {
        raw_verts.resize(total_raw_vertices);
    } catch (const std::bad_alloc&) {
        throw std::runtime_error("OOM: Failed to allocate memory for STL geometry.");
    }

    #pragma pack(push, 1)
    struct StlTriangle { float n[3]; float v1[3]; float v2[3]; float v3[3]; uint16_t attr; };
    #pragma pack(pop)

    const size_t CHUNK_SIZE = 5000;
    std::vector<StlTriangle> buffer(CHUNK_SIZE);

    size_t triangles_read = 0;
    uint32_t v_idx = 0;

    // 1. Read all vertices sequentially (Maximum disk/memory bandwidth)
    while (triangles_read < num_triangles) {
        size_t batch_size = std::min(CHUNK_SIZE, static_cast<size_t>(num_triangles - triangles_read));
        file.read(reinterpret_cast<char*>(buffer.data()), batch_size * 50);

        if (!file && file.gcount() < static_cast<std::streamsize>(batch_size * 50)) {
            throw std::runtime_error("Unexpected EOF: STL file is corrupted or truncated.");
        }

        for (size_t i = 0; i < batch_size; ++i) {
            // Note: "+ 0.0f" normalizes IEEE-754 -0.0f to 0.0f silently in 1 cycle
            raw_verts[v_idx] = {buffer[i].v1[0] + 0.0f, buffer[i].v1[1] + 0.0f, buffer[i].v1[2] + 0.0f, v_idx}; v_idx++;
            raw_verts[v_idx] = {buffer[i].v2[0] + 0.0f, buffer[i].v2[1] + 0.0f, buffer[i].v2[2] + 0.0f, v_idx}; v_idx++;
            raw_verts[v_idx] = {buffer[i].v3[0] + 0.0f, buffer[i].v3[1] + 0.0f, buffer[i].v3[2] + 0.0f, v_idx}; v_idx++;
        }
        triangles_read += batch_size;
    }

    // 2. Sort in contiguous memory (Blisteringly fast, highly cache localized)
    std::sort(raw_verts.begin(), raw_verts.end());

    // 3. Sweep and Remap
    std::vector<uint32_t> remap(total_raw_vertices);
    std::vector<Eigen::Vector3f> unique_vertices;
    unique_vertices.reserve(num_triangles / 2); // educated guess

    unique_vertices.emplace_back(raw_verts[0].x, raw_verts[0].y, raw_verts[0].z);
    remap[raw_verts[0].original_id] = 0;
    
    uint32_t current_unique_id = 0;
    for (size_t i = 1; i < total_raw_vertices; ++i) {
        if (!(raw_verts[i] == raw_verts[i - 1])) {
            current_unique_id++;
            unique_vertices.emplace_back(raw_verts[i].x, raw_verts[i].y, raw_verts[i].z);
        }
        remap[raw_verts[i].original_id] = current_unique_id;
    }

    // 4. Construct Final FastMesh
    FastMesh mesh;
    mesh.vertices.resize(unique_vertices.size(), 3);
    mesh.triangles.resize(num_triangles, 3);

    #pragma omp parallel for
    for (int i = 0; i < (int)unique_vertices.size(); ++i) {
        mesh.vertices.row(i) = unique_vertices[i];
    }

    #pragma omp parallel for
    for (int i = 0; i < (int)num_triangles; ++i) {
        mesh.triangles.row(i) << remap[i * 3 + 0], remap[i * 3 + 1], remap[i * 3 + 2];
    }

    #ifdef MULTISPHERE_DEBUG
        std::cout << "[IO] Loaded STL (Sort-and-Sweep). Vertices welded: " 
                  << total_raw_vertices << " -> " << unique_vertices.size() << std::endl;
    #endif

    return mesh;
}

/**
 * @brief Exports a SpherePack to CSV.
 * @param sp SpherePack to export.
 * @param path Output CSV file path.
 */
inline void export_to_csv(const SpherePack& sp, const std::string& path) {
    std::ofstream file(path);
    // Add a buffer to speed up I/O
    std::vector<char> buffer(65536); 
    file.rdbuf()->pubsetbuf(buffer.data(), buffer.size());

    file << "x,y,z,radius\n";
    for (int i = 0; i < sp.num_spheres(); ++i) {
        file << sp.centers(i, 0) << "," 
             << sp.centers(i, 1) << "," 
             << sp.centers(i, 2) << "," 
             << sp.radii(i) << "\n";
    }
}

/**
 * @brief Exports a SpherePack to legacy VTK format.
 * @param sp SpherePack to export.
 * @param path Output VTK file path.
 */
inline void export_to_vtk(const SpherePack& sp, const std::string& path) {
    std::ofstream f(path);
    size_t n = sp.num_spheres();

    // Legacy VTK format is the most efficient for point clouds
    f << "# vtk DataFile Version 3.0\n"
      << "SpherePack Centers\nASCII\nDATASET POLYDATA\n"
      << "POINTS " << n << " float\n";

    for (int i = 0; i < n; ++i) {
        f << (float)sp.centers(i, 0) << " " 
          << (float)sp.centers(i, 1) << " " 
          << (float)sp.centers(i, 2) << "\n";
    }

    f << "VERTICES " << n << " " << n * 2 << "\n";
    for (int i = 0; i < n; ++i) f << "1 " << i << "\n";

    f << "POINT_DATA " << n << "\nSCALARS radius float 1\nLOOKUP_TABLE default\n";
    for (int i = 0; i < n; ++i) f << (float)sp.radii(i) << "\n";
}

/**
 * @brief Exports a VoxelGrid distance transform to a Legacy VTK file.
 * Now templatized to handle VoxelGrid<T> regardless of the underlying data type.
 * @tparam T VoxelGrid data type.
 * @param grid VoxelGrid to export.
 * @param path Output VTK file path.
 */
template <typename T>
inline void export_voxel_grid_to_vtk(const VoxelGrid<T>& grid, const std::string& path) {
    std::ofstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file for writing: " + path);
    }

    int nx = grid.shape[0];
    int ny = grid.shape[1];
    int nz = grid.shape[2];

    file << "# vtk DataFile Version 3.0\n";
    file << "VoxelGrid Distance Transform\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";

    file << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";
    file << "ORIGIN " << grid.origin[0] << " " << grid.origin[1] << " " << grid.origin[2] << "\n";
    file << "SPACING " << grid.voxel_size << " " << grid.voxel_size << " " << grid.voxel_size << "\n";

    file << "POINT_DATA " << (nx * ny * nz) << "\n";
    
    std::string type_name = "float";
    if constexpr (std::is_same_v<T, double>) type_name = "double";
    if constexpr (std::is_same_v<T, int>)    type_name = "int";

    file << "SCALARS distance_transform " << type_name << " 1\n";
    file << "LOOKUP_TABLE default\n";

    for (const auto& value : grid.data) {
        file << value << "\n";
    }

    file.close();
}

/**
 * @brief Loads a VoxelGrid<bool> from a .npy file.
 * @param path Path to .npy file.
 * @param voxel_size Physical voxel size.
 * @return VoxelGrid<uint8_t> loaded from file.
 */
inline VoxelGrid<uint8_t> load_voxels_from_npy(const std::string& path, double voxel_size) {
    cnpy::NpyArray arr = cnpy::npy_load(path);
    if (arr.shape.size() != 3) throw std::runtime_error("Array must be 3D");

    VoxelGrid<uint8_t> grid(arr.shape[0], arr.shape[1], arr.shape[2], voxel_size);
    
    // Treat numpy boolean data natively as 1-byte uint8_t
    if (arr.word_size == 1) {
        uint8_t* loaded_data = arr.data<uint8_t>();
        std::copy(loaded_data, loaded_data + grid.data.size(), grid.data.begin());
    } else {
        throw std::runtime_error("Unsupported NPY word size for binary grid loading.");
    }
    
    return grid;
}

/**
 * @brief Saves a VoxelGrid to a .npy file.
 * @tparam T VoxelGrid data type.
 * @param path Output .npy file path.
 * @param grid VoxelGrid to save.
 */
template<typename T>
inline void save_voxels_to_npy(const std::string& path, const VoxelGrid<T>& grid) {
    std::vector<size_t> shape = {
        static_cast<size_t>(grid.shape[0]), 
        static_cast<size_t>(grid.shape[1]), 
        static_cast<size_t>(grid.shape[2])
    };
    cnpy::npy_save(path, grid.data.data(), shape, "w");
}

/**
 * @brief Saves a FastMesh to a binary STL file.
 * @param mesh FastMesh to save.
 * @param output_path Output STL file path.
 */
inline void save_mesh_to_stl(const FastMesh& mesh, const std::string& output_path) {
    std::ofstream file(output_path, std::ios::binary);
    if (!file) throw std::runtime_error("Could not open file for writing: " + output_path);

    char header[80] = {0};
    std::string title = "Debug Export";
    std::copy(title.begin(), title.end(), header);
    file.write(header, 80);

    uint32_t n_tris = static_cast<uint32_t>(mesh.triangles.rows());
    file.write(reinterpret_cast<const char*>(&n_tris), 4);

    for (uint32_t i = 0; i < n_tris; ++i) {
        float normal[3] = {0.0f, 0.0f, 0.0f};
        file.write(reinterpret_cast<char*>(normal), 12);

        for (int v = 0; v < 3; ++v) {
            int v_idx = mesh.triangles(i, v);
            float x = mesh.vertices(v_idx, 0);
            float y = mesh.vertices(v_idx, 1);
            float z = mesh.vertices(v_idx, 2);
            float vert_data[3] = {x, y, z};
            file.write(reinterpret_cast<char*>(vert_data), 12);
        }

        uint16_t attr = 0;
        file.write(reinterpret_cast<char*>(&attr), 2);
    }
    #ifdef MULTISPHERE_DEBUG
        std::cout << "[IO] Saved debug mesh to: " << output_path << std::endl;
    #endif

}


} // namespace MSS

#endif // MULTISPHERE_IO_HPP
