#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <fstream>
#include <random>
#include <cmath>
#include <Eigen/Dense>

#include "GEMSS/GEMSS-interface.h"
#include "GEMSS/GEMSS_split.hpp"

using namespace GEMSS;
namespace fs = std::filesystem;

// Computes the geometric intersection volume between two SpherePacks
float compute_fragment_overlap(const SpherePack& frag1, const SpherePack& frag2) {
    float total_overlap_volume = 0.0f;
    
    for (int i = 0; i < frag1.centers.rows(); ++i) {
        Eigen::Vector3f c1 = frag1.centers.row(i).transpose();
        float r1 = frag1.radii(i);
        
        for (int j = 0; j < frag2.centers.rows(); ++j) {
            Eigen::Vector3f c2 = frag2.centers.row(j).transpose();
            float r2 = frag2.radii(j);
            
            float d = (c1 - c2).norm();
            
            if (d >= r1 + r2) {
                continue; // No overlap
            } else if (d <= std::abs(r1 - r2)) {
                // One completely inside the other
                float min_r = std::min(r1, r2);
                total_overlap_volume += (4.0f / 3.0f) * M_PI * std::pow(min_r, 3);
            } else {
                // Partial overlap
                float v = (M_PI * std::pow(r1 + r2 - d, 2) * (std::pow(d, 2) + 2 * d * (r1 + r2) - 3 * std::pow(r1 - r2, 2))) / (12.0f * d);
                total_overlap_volume += v;
            }
        }
    }
    return total_overlap_volume;
}

// Generates a random unit normal vector
Eigen::Vector3f generate_random_normal(std::mt19937& gen) {
    std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
    Eigen::Vector3f normal;
    do {
        normal = Eigen::Vector3f(dist(gen), dist(gen), dist(gen));
    } while (normal.squaredNorm() < 0.01f);
    return normal.normalized();
}

// Generates a random plane point guaranteed to intersect the SpherePack
Eigen::Vector3f generate_intersecting_point(const SpherePack& sp, const Eigen::Vector3f& normal, std::mt19937& gen) {
    float min_proj = std::numeric_limits<float>::max();
    float max_proj = std::numeric_limits<float>::lowest();

    for (int i = 0; i < sp.centers.rows(); ++i) {
        float proj = normal.dot(sp.centers.row(i).transpose());
        float r = sp.radii(i);
        if (proj - r < min_proj) min_proj = proj - r;
        if (proj + r > max_proj) max_proj = proj + r;
    }

    // Shrink bounds slightly to avoid glancing cuts at the absolute extremeties
    float buffer = (max_proj - min_proj) * 0.1f; 
    std::uniform_real_distribution<float> dist(min_proj + buffer, max_proj - buffer);
    
    return normal * dist(gen);
}

int main(int argc, char** argv) {
    std::string dataset_dir = "./samples"; // Directory containing STL files
    
    std::random_device rd;
    std::mt19937 gen(rd());

    // Define tuning ranges for the study
    std::vector<float> parent_offsets = {0.0f, 0.25f, 0.5f, 0.75f, 1.0f};
    std::vector<float> plane_offsets = {0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f};
    std::vector<float> divs = {200,250};

for (int div :divs){
    GEMSS::MultisphereConfig base_config;
    base_config.div = div;
    base_config.compute_physics = 1;
    base_config.conserve_mass = false; // CRITICAL: Disabled to measure natural mass loss
    base_config.show_progress = false;

    std::string output_csv = std::to_string(div) + "tuning_results.csv";

    std::ofstream out_file(output_csv);
    out_file << "filename,parent_rep_offset,fracture_plane_offset,parent_mass,fragments_mass,mass_ratio,overlap_volume\n";

    


    for (const auto& entry : fs::directory_iterator(dataset_dir)) {
        if (entry.path().extension() != ".stl") continue;
        
        std::string filename = entry.path().filename().string();
        STLMesh mesh = load_mesh(entry.path().string());
        if (mesh.is_empty()) continue;
        SpherePack parent_sp = multisphere_from_mesh(mesh, base_config);
        
        std::cout<< "Starting " << div << " on " << filename <<std::endl; 

        if (parent_sp.mass <= 0.0f) {
            std::cout<< filename << " sp failed!"<<std::endl;
            continue;
        }

        Eigen::Vector3f normal = generate_random_normal(gen);
        Eigen::Vector3f point = generate_intersecting_point(parent_sp, normal, gen);

        for (float p_offset : parent_offsets) {
            for (float f_offset : plane_offsets) {
                
                GEMSS::MultisphereConfig study_config = base_config;
                study_config.parent_representation_offset = p_offset;
                study_config.fracture_plane_offset = f_offset;

                auto [fragments, labeled_grid, area] = split_and_compute_surface_sp(parent_sp, normal, point, study_config);

                float total_frag_mass = 0.0f;
                for (const auto& frag : fragments) {
                    total_frag_mass += frag.mass;
                }

                float mass_ratio = total_frag_mass / parent_sp.mass;

                float overlap_vol = 0.0f;
                if (fragments.size() >= 2) {
                    overlap_vol = compute_fragment_overlap(fragments[0], fragments[1]);
                }

                out_file << filename << ","
                         << p_offset << ","
                         << f_offset << ","
                         << parent_sp.mass << ","
                         << total_frag_mass << ","
                         << mass_ratio << ","
                         << overlap_vol << "\n";
            
            
                //std::cout<<filename<<"_"<<p_offset<<"_"<<f_offset<<std::endl;
            }
        }
    }

    out_file.close();
}
    return 0;
}