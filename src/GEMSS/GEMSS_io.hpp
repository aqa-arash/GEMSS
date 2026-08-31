/**
 * @file multisphere_io.hpp
 * @brief Input/output utilities for multisphere-cpp library.
 *
 * Provides mesh loading, sphere pack export, voxel grid export, and STL/NPY file utilities.
 *
 * @author Arash Moradian
 * @date 2026-03-09
 */

#ifndef GEMSS_IO_HPP
#define GEMSS_IO_HPP

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <filesystem>
#include <system_error>
#include <Eigen/Core>
#include "GEMSS_datatypes.hpp"

namespace GEMSS {

/**
 * @brief Loads a mesh from a binary STL file.
 * @param path Path to STL file.
 * @return STLMesh structure.
 */
inline STLMesh load_mesh(const std::string& path) {
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
    if (ec || physical_size < 84 + (static_cast<uintmax_t>(num_triangles) * 50)) {
        throw std::runtime_error("Malformed Binary STL: Cannot read file or file too small for claimed triangle count.");
    }

    if (num_triangles == 0) return STLMesh();

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

    // 2. Sort in contiguous memory for vertex welding (Exploits spatial locality and SIMD during sorting)
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

    // 4. Construct Final STLMesh
    STLMesh mesh;
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

    // VTK format
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
 * @brief Saves a STLMesh to a binary STL file.
 * @param mesh STLMesh to save.
 * @param output_path Output STL file path.
 */
inline void save_mesh_to_stl(const STLMesh& mesh, const std::string& output_path) {
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


/**
 * @brief Print detailed information about a SpherePack to the console.
 */
inline void print_sphere_pack_info(const SpherePack& sp) {
    std::cout << "Sphere Pack Info:" << std::endl;
    std::cout << "      Number of spheres: " << sp.num_spheres() << std::endl;
    std::cout << "      Max radius: " << sp.max_radius() << " units" << std::endl;
    std::cout << "      Min radius: " << sp.min_radius() << " units" << std::endl;
    std::cout << "      Accuracy: " << sp.precision << " units^3" << std::endl;
    std::cout << "      Mass of union: " << sp.mass << " units^3" << std::endl;
    std::cout << "      Center of mass: " << sp.center_of_mass.transpose() << " units" << std::endl;
    std::cout << "      Bounding radius: " << sp.bounding_radius << " units" << std::endl;
    std::cout << "      Principal moments: " << sp.principal_moments.transpose() << " units^5" << std::endl;
    std::cout << "      Principal axes:\n" << sp.principal_axes << std::endl;

}

/**
 * @brief Export sphere pack to JSON file 
 */
inline void export_sphere_pack_json(const SpherePack& sp, const std::string& output_path) {
    // check if output has path and make it 
    std::filesystem::path outp(output_path);
    if (!outp.parent_path().empty()) std::filesystem::create_directories(outp.parent_path());

    //open the file 
    std::ofstream f(output_path, std::ios::out | std::ios::trunc);
    if (!f.is_open()) throw std::runtime_error("Could not open file for writing: " + output_path);

    // helper to maintain correct json format 
    auto write_num = [&](double v) {
        if (std::isfinite(v)) f << v; else f << "null";
    };


    f << std::fixed << std::setprecision(9); //for float, if needed switch to double = 17
    f << "{\n";

    f << "  \"n_spheres\": " << sp.num_spheres() << ",\n";

    // centers
    f << "  \"centers\": [\n";
    for (int i = 0; i < (int)sp.centers.rows(); ++i) {
        f << "    [";
        write_num(sp.centers(i,0)); f << ", ";
        write_num(sp.centers(i,1)); f << ", ";
        write_num(sp.centers(i,2));
        f << "]";
        if (i + 1 < sp.centers.rows()) f << ",\n"; else f << "\n";
    }
    f << "  ],\n";

    // radii
    f << "  \"radii\": [";
    for (int i = 0; i < (int)sp.radii.size(); ++i) {
        write_num(sp.radii(i));
        if (i + 1 < sp.radii.size()) f << ", ";
    }
    f << "],\n";

    // simple fields
    f << "  \"is_2d\": " << (sp.is_2d ? "true" : "false") << ",\n";
    f << "  \"precision\": "; write_num(sp.precision); f << ",\n";
    f << "  \"density\": "; write_num(sp.density); f << ",\n";
    f << "  \"mass\": "; write_num(sp.mass); f << ",\n";
    f << "  \"bounding_radius\": "; write_num(sp.bounding_radius); f << ",\n";

    // center_of_mass
    f << "  \"center_of_mass\": [";
    write_num(sp.center_of_mass.x()); f << ", ";
    write_num(sp.center_of_mass.y()); f << ", ";
    write_num(sp.center_of_mass.z()); f << "],\n";

    // inertia_tensor (row-major)
    f << "  \"inertia_tensor\": [";
    for (int r = 0; r < 3; ++r) for (int c = 0; c < 3; ++c) {
        write_num(sp.inertia_tensor(r,c));
        if (!(r==2 && c==2)) f << ", ";
    }
    f << "],\n";

    // principal_axes (row-major)
    f << "  \"principal_axes\": [";
    for (int r = 0; r < 3; ++r) for (int c = 0; c < 3; ++c) {
        write_num(sp.principal_axes(r,c));
        if (!(r==2 && c==2)) f << ", ";
    }
    f << "],\n";

    // principal_moments
    f << "  \"principal_moments\": [";
    write_num(sp.principal_moments.x()); f << ", ";
    write_num(sp.principal_moments.y()); f << ", ";
    write_num(sp.principal_moments.z());
    f << "]\n";

    f << "}\n";
    f.close();
}


/**
 * @brief imports spherepack from , csv (center and radii only)
*/
inline SpherePack import_csv(const std::string& path) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("File not found: " + path);

    std::string line;
    std::getline(file, line); 

    std::vector<Eigen::Vector3f> centers;
    std::vector<float> radii;

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string token;
        float x, y, z, r;

        std::getline(ss, token, ','); x = std::stof(token);
        std::getline(ss, token, ','); y = std::stof(token);
        std::getline(ss, token, ','); z = std::stof(token);
        std::getline(ss, token, ','); r = std::stof(token);

        centers.emplace_back(x, y, z);
        radii.push_back(r);
    }

    SpherePack sp(Eigen::MatrixX3f(centers.size(), 3), Eigen::VectorXf(radii.size()));
    for (size_t i = 0; i < centers.size(); ++i) {
        sp.centers.row(i) = centers[i];
        sp.radii(i) = radii[i];
    }

    std::cout << "Import made from CSV format." << std::endl;
    return sp;
}

/**
 * @brief import center and radii from csv 
 */
inline SpherePack import_vtk(const std::string& path) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("File not found: " + path);

    std::string token;
    size_t num_points = 0;

    while (file >> token) {
        if (token == "POINTS") {
            file >> num_points >> token; 
            break;
        }
    }

    SpherePack sp(Eigen::MatrixX3f(num_points, 3), Eigen::VectorXf(num_points));
    for (size_t i = 0; i < num_points; ++i) {
        file >> sp.centers(i, 0) >> sp.centers(i, 1) >> sp.centers(i, 2);
    }

    while (file >> token) {
        if (token == "SCALARS") {
            file >> token >> token >> token; 
            file >> token >> token; 
            break;
        }
    }

    for (size_t i = 0; i < num_points; ++i) file >> sp.radii(i);

    std::cout << "Import made from VTK format." << std::endl;
    return sp;
}

/**
 * @brief Hardcoded Json parser to remove library dependency 
 */

 inline SpherePack import_json(const std::string& path) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("File not found: " + path);
    
    // Load entire file into memory for fast string manipulation
    std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

    // Extract sequential numeric values associated with a specific JSON key
    auto extract_numbers = [&](const std::string& key, size_t count) -> std::vector<float> {
        std::vector<float> vals;
        size_t pos = content.find("\"" + key + "\"");
        if (pos == std::string::npos) {
            std::cout << "[Warning] '" << key << "' missing.\n";
            return vals;
        }
        
        pos = content.find(':', pos);
        if (pos == std::string::npos) return vals;
        
        pos++; // Advance index past the colon separator

        std::string token;
        while (pos < content.size() && vals.size() < count) {
            char c = content[pos];
            if (c == '"') break; 
            
            if (std::isspace(c) || c == ',' || c == '[' || c == ']') {
                if (!token.empty()) {
                    if (token == "null") vals.push_back(std::numeric_limits<float>::quiet_NaN());
                    else vals.push_back(std::stof(token));
                    token.clear();
                }
            } else {
                token += c;
            }
            pos++;
        }
        
        if (!token.empty() && vals.size() < count) {
            if (token == "null") vals.push_back(std::numeric_limits<float>::quiet_NaN());
            else vals.push_back(std::stof(token));
        }
        return vals;
    };

    // Locate key and determine boolean state
    auto extract_bool = [&](const std::string& key) -> bool {
        size_t pos = content.find("\"" + key + "\"");
        if (pos == std::string::npos) return false;
        pos = content.find(':', pos);
        size_t t = content.find("true", pos);
        size_t f = content.find("false", pos);
        return (t < f); 
    };

    SpherePack sp;

    std::vector<float> n_val = extract_numbers("n_spheres", 1);
    if (n_val.empty()) return sp;
    size_t n = static_cast<size_t>(n_val[0]);

    std::vector<float> c_vals = extract_numbers("centers", n * 3);
    std::vector<float> r_vals = extract_numbers("radii", n);

    if (c_vals.size() == n * 3 && r_vals.size() == n) {
        sp.centers.resize(n, 3);
        sp.radii.resize(n);
        for (size_t i = 0; i < n; ++i) {
            sp.centers(i, 0) = c_vals[i * 3 + 0];
            sp.centers(i, 1) = c_vals[i * 3 + 1];
            sp.centers(i, 2) = c_vals[i * 3 + 2];
            sp.radii(i) = r_vals[i];
        }
    }

    sp.is_2d = extract_bool("is_2d");

    auto set_scalar = [&](const std::string& key, float& target) {
        std::vector<float> v = extract_numbers(key, 1);
        if (!v.empty()) target = v[0];
    };

    set_scalar("precision", sp.precision);
    set_scalar("density", sp.density);
    set_scalar("mass", sp.mass);
    set_scalar("bounding_radius", sp.bounding_radius);

    auto set_vec = [&](const std::string& key, Eigen::Vector3f& target) {
        std::vector<float> v = extract_numbers(key, 3);
        if (v.size() == 3) target << v[0], v[1], v[2];
    };

    set_vec("center_of_mass", sp.center_of_mass);
    set_vec("principal_moments", sp.principal_moments);

    auto set_mat = [&](const std::string& key, Eigen::Matrix3f& target) {
        std::vector<float> v = extract_numbers(key, 9);
        if (v.size() == 9) {
            for (int r = 0; r < 3; ++r) 
                for (int c = 0; c < 3; ++c) 
                    target(r, c) = v[r * 3 + c];
        }
    };

    set_mat("inertia_tensor", sp.inertia_tensor);
    set_mat("principal_axes", sp.principal_axes);

    std::cout << "Import made from JSON format." << std::endl;
    return sp;
}


/**
 * @brief wrapper for the import functions 
 */
inline SpherePack import_sphere_pack(const std::string& path) {
    std::filesystem::path file_path(path);
    std::string ext = file_path.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    if (ext == ".csv") return import_csv(path);
    if (ext == ".vtk") return import_vtk(path);
    if (ext == ".json") return import_json(path);

    throw std::runtime_error("Unsupported file extension: " + ext);
}
} // namespace MSS

#endif // GEMSS_IO_HPP