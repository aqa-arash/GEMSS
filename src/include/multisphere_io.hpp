#ifndef MULTISPHERE_IO_HPP
#define MULTISPHERE_IO_HPP

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <Eigen/Dense>
#include "include/cnpy.h"
#include "multisphere_datatypes.hpp"

inline FastMesh load_mesh_fast(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) throw std::runtime_error("File not found: " + path);

    // Skip Header
    file.seekg(80);
    uint32_t num_triangles;
    file.read(reinterpret_cast<char*>(&num_triangles), 4);

    std::vector<Eigen::Vector3f> unique_vertices;
    std::vector<Eigen::Vector3i> face_indices;
    
    // Reserve memory to prevent reallocations (Critical for speed)
    unique_vertices.reserve(num_triangles / 2); 
    face_indices.reserve(num_triangles);

    // Use Unordered Map for O(1) lookups
    std::unordered_map<VertexKey, int, VertexKeyHash> vertex_map;
    vertex_map.reserve(num_triangles * 2);

    const size_t CHUNK_SIZE = 5000;
    
    #pragma pack(push, 1)
    struct StlTriangle { float n[3]; float v1[3]; float v2[3]; float v3[3]; uint16_t attr; };
    #pragma pack(pop)
    
    std::vector<StlTriangle> buffer(CHUNK_SIZE);

    size_t triangles_read = 0;
    while (triangles_read < num_triangles) {
        size_t batch_size = std::min(CHUNK_SIZE, (size_t)(num_triangles - triangles_read));
        file.read(reinterpret_cast<char*>(buffer.data()), batch_size * 50);

        for (size_t i = 0; i < batch_size; ++i) {
            int indices[3];
            float* raw_verts[3] = { buffer[i].v1, buffer[i].v2, buffer[i].v3 };

            for (int k = 0; k < 3; ++k) {
                VertexKey key = { raw_verts[k][0], raw_verts[k][1], raw_verts[k][2] };

                // O(1) Lookup
                auto it = vertex_map.find(key);
                if (it != vertex_map.end()) {
                    indices[k] = it->second;
                } else {
                    int new_idx = static_cast<int>(unique_vertices.size());
                    unique_vertices.emplace_back(key.x, key.y, key.z);
                    vertex_map.insert({key, new_idx}); // slightly faster than operator[]
                    indices[k] = new_idx;
                }
            }
            face_indices.emplace_back(indices[0], indices[1], indices[2]);
        }
        triangles_read += batch_size;
    }

    FastMesh mesh;
    // Map directly to Eigen (no copy loop needed usually, but safe to keep explicit)
    mesh.vertices.resize(unique_vertices.size(), 3);
    mesh.triangles.resize(face_indices.size(), 3);

    // Parallel copy if mesh is huge
    #pragma omp parallel for
    for (int i = 0; i < (int)unique_vertices.size(); ++i) mesh.vertices.row(i) = unique_vertices[i];
    #pragma omp parallel for
    for (int i = 0; i < (int)face_indices.size(); ++i) mesh.triangles.row(i) = face_indices[i];

    std::cout << "[IO] Loaded STL. Vertices welded: " << (num_triangles*3) << " -> " << unique_vertices.size() << std::endl;
    return mesh;
}


void export_to_csv(const SpherePack& sp, const std::string& path) {
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

void export_to_vtk(const SpherePack& sp, const std::string& path) {
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
 */
template <typename T>
void export_voxel_grid_to_vtk(const VoxelGrid<T>& grid, const std::string& path) {
    std::ofstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file for writing: " + path);
    }

    // Accessing members (Ensure these are public or have getters in your VoxelGrid class)
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
    
    // Logic: Identify VTK type based on T
    std::string type_name = "float";
    if constexpr (std::is_same_v<T, double>) type_name = "double";
    if constexpr (std::is_same_v<T, int>)    type_name = "int";

    file << "SCALARS distance_transform " << type_name << " 1\n";
    file << "LOOKUP_TABLE default\n";

    // Flattened data output
    // Assumption: grid.data is a std::vector<T> or similar container
    for (const auto& value : grid.data) {
        file << value << "\n";
    }

    file.close();
}


VoxelGrid<bool> load_voxels_from_npy(const std::string& path, double voxel_size) {
    cnpy::NpyArray arr = cnpy::npy_load(path);
    if (arr.shape.size() != 3) throw std::runtime_error("Array must be 3D");

    VoxelGrid<bool> grid(arr.shape[0], arr.shape[1], arr.shape[2], voxel_size);
    
    // Direct memory copy (memcpy) is the fastest way to load the data
    bool* loaded_data = arr.data<bool>();
    std::copy(loaded_data, loaded_data + grid.data.size(), grid.data.begin());
    
    return grid;
}


template<typename T>
void save_voxels_to_npy(const std::string& path, const VoxelGrid<T>& grid) {
    // 1. Define the shape of the 3D array for NumPy (Z, Y, X)
    // Note: NumPy usually expects C-style (Row-Major) layout
    std::vector<size_t> shape = {
        static_cast<size_t>(grid.shape[0]), 
        static_cast<size_t>(grid.shape[1]), 
        static_cast<size_t>(grid.shape[2])
    };

    // 2. Use npy_save
    // Arguments: Filename, Data Pointer, Shape Vector, Mode ("w" overrides, "a" appends)
    cnpy::npy_save(path, grid.data.data(), shape, "w");
}


inline void save_mesh_to_stl(const FastMesh& mesh, const std::string& output_path) {
    std::ofstream file(output_path, std::ios::binary);
    if (!file) throw std::runtime_error("Could not open file for writing: " + output_path);

    // 1. Write Header (80 bytes)
    char header[80] = {0};
    std::string title = "Debug Export";
    std::copy(title.begin(), title.end(), header);
    file.write(header, 80);

    // 2. Write Triangle Count (4 bytes)
    uint32_t n_tris = static_cast<uint32_t>(mesh.triangles.rows());
    file.write(reinterpret_cast<const char*>(&n_tris), 4);

    // 3. Write Triangles (50 bytes each)
    // Structure: Normal(12) + V1(12) + V2(12) + V3(12) + Attr(2)
    for (uint32_t i = 0; i < n_tris; ++i) {
        // A. Normal - We just write zero (0,0,0) as we don't store normals in FastMesh
        // Most viewers (MeshLab, Paraview) will auto-calculate them.
        float normal[3] = {0.0f, 0.0f, 0.0f};
        file.write(reinterpret_cast<char*>(normal), 12);

        // B. Vertices
        for (int v = 0; v < 3; ++v) {
            int v_idx = mesh.triangles(i, v);
            // Access Eigen row
            float x = mesh.vertices(v_idx, 0);
            float y = mesh.vertices(v_idx, 1);
            float z = mesh.vertices(v_idx, 2);
            
            float vert_data[3] = {x, y, z};
            file.write(reinterpret_cast<char*>(vert_data), 12);
        }

        // C. Attribute byte count (2 bytes)
        uint16_t attr = 0;
        file.write(reinterpret_cast<char*>(&attr), 2);
    }
    
    std::cout << "[IO] Saved debug mesh to: " << output_path << std::endl;
}




#endif