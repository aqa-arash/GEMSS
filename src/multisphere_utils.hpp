#ifndef MULTISPHERE_UTILS_HPP
#define MULTISPHERE_UTILS_HPP

#ifdef HAVE_MANIFOLD

#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include "manifold/manifold.h"
#include "nanoflann.hpp"
#include "multisphere_datatypes.hpp"

using namespace manifold;
namespace utils {
    // --- Helper: Mesh to Manifold ---
    inline Manifold mesh_to_manifold(const std::vector<Eigen::Vector3f>& vertices, 
                                     const std::vector<std::array<int, 3>>& faces) {
        MeshGL mesh_gl;
        mesh_gl.numProp = 3;
        for (const auto& v : vertices) {
            mesh_gl.vertPos.push_back(v.x());
            mesh_gl.vertPos.push_back(v.y());
            mesh_gl.vertPos.push_back(v.z());
        }
        for (const auto& f : faces) {
            mesh_gl.triVerts.push_back(f[0]);
            mesh_gl.triVerts.push_back(f[1]);
            mesh_gl.triVerts.push_back(f[2]);
        }
        return Manifold(mesh_gl);
    }

    // --- Dice Coefficient ---
    inline double compute_dice_coefficient(const Manifold& m1, const Manifold& m2) {
        double v1 = m1.GetProperties().volume;
        double v2 = m2.GetProperties().volume;
        if (v1 + v2 <= 0) return 0.0;

        Manifold intersection = m1 ^ m2; 
        return (2.0 * intersection.GetProperties().volume) / (v1 + v2) * 100.0;
    }

    // --- Create Multisphere Mesh ---
    inline Manifold create_multisphere_mesh(const SpherePack& sp, int resolution = 4) {
        if (sp.num_spheres() == 0) throw std::runtime_error("No spheres provided.");
        std::vector<Manifold> sphere_list;
        for (size_t i = 0; i < sp.num_spheres(); ++i) {
            Manifold s = Manifold::Sphere(sp.radii(i), resolution * 4);
            s = s.Translate({(float)sp.centers(i,0), (float)sp.centers(i,1), (float)sp.centers(i,2)});
            sphere_list.push_back(s);
        }
        return Manifold::BatchBoolean(sphere_list, OpType::Add);
    }

    // --- Adjust Spheres to Boundary ---
    inline SpherePack adjust_spheres_to_stl_boundary(const SpherePack& sp, 
                                                     const std::vector<Eigen::Vector3f>& mesh_vertices) {
        using my_kd_tree_t = nanoflann::KDTreeEigenMatrixAdaptor<Eigen::Matrix<float, Eigen::Dynamic, 3>>;
        Eigen::Matrix<float, Eigen::Dynamic, 3> mat(mesh_vertices.size(), 3);
        for(size_t i=0; i<mesh_vertices.size(); ++i) mat.row(i) = mesh_vertices[i];

        my_kd_tree_t index(3, std::cref(mat), 10);
        SpherePack adjusted = sp;

        for (size_t i = 0; i < adjusted.num_spheres(); ++i) {
            float query_pt[3] = {(float)adjusted.centers(i,0), (float)adjusted.centers(i,1), (float)adjusted.centers(i,2)};
            size_t ret_index;
            float out_dist_sq;
            nanoflann::KNNResultSet<float> resultSet(1);
            resultSet.init(&ret_index, &out_dist_sq);
            index.index->findNeighbors(resultSet, &query_pt[0]);

            float distance = std::sqrt(out_dist_sq);
            if (distance < adjusted.radii(i)) adjusted.radii(i) = distance;
        }
        return adjusted;
    }
}


namespace fs = std::filesystem;

/**
 * Export a multisphere SpherePack as an STL mesh file.
 * Logic: Boolean Union (via Manifold) -> Binary STL Export.
 */
void export_sphere_pack_to_stl(
    const SpherePack& sphere_pack,
    const std::string& path,
    int resolution = 4
) {
    // 1. Extension and Directory Validation
    fs::path export_path(path);
    if (export_path.extension() != ".stl" && export_path.extension() != ".STL") {
        throw std::invalid_argument("Output path must end with '.stl'.");
    }

    fs::path directory = export_path.parent_path();
    if (!directory.empty() && !fs::exists(directory)) {
        throw std::runtime_error("Directory for STL export does not exist: " + directory.string());
    }

    // 2. Generate the Union Mesh (Only if Manifold is enabled)
#ifdef HAVE_MANIFOLD
    // create_multisphere_mesh returns a manifold::Manifold object
    manifold::Manifold result_mesh = utils::create_multisphere_mesh(sphere_pack, resolution);
    
    // Get MeshGL structure (standard manifold output format)
    manifold::MeshGL mesh_out = result_mesh.GetMeshGL();

    // 3. Binary STL Export
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs) throw std::runtime_error("Could not open file for writing: " + path);

    // Write 80-byte header
    char header[80] = "Multisphere generated STL";
    ofs.write(header, 80);

    // Write number of triangles
    uint32_t num_triangles = static_cast<uint32_t>(mesh_out.triVerts.size() / 3);
    ofs.write(reinterpret_cast<const char*>(&num_triangles), 4);

    // Each triangle: 3*4 (normal) + 3*3*4 (vertices) + 2 (attr) = 50 bytes
    for (size_t i = 0; i < mesh_out.triVerts.size(); i += 3) {
        float dummy_normal[3] = {0.0f, 0.0f, 0.0f}; // STL normals can be zero
        ofs.write(reinterpret_cast<const char*>(dummy_normal), 12);

        for (int j = 0; j < 3; ++j) {
            uint32_t vert_idx = mesh_out.triVerts[i + j];
            float v[3] = {
                mesh_out.vertPos[vert_idx * 3],
                mesh_out.vertPos[vert_idx * 3 + 1],
                mesh_out.vertPos[vert_idx * 3 + 2]
            };
            ofs.write(reinterpret_cast<const char*>(v), 12);
        }

        uint16_t attr = 0;
        ofs.write(reinterpret_cast<const char*>(&attr), 2);
    }
#else
    throw std::runtime_error("export_sphere_pack_to_stl requires Manifold library.");
#endif
}

#endif // HAVE_MANIFOLD
#endif // MULTISPHERE_UTILS_HPP