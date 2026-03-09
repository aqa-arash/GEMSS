// multisphere_reconstruction_helpers.hpp
#ifndef MULTISPHERE_RECONSTRUCTION_HELPERS_HPP
#define MULTISPHERE_RECONSTRUCTION_HELPERS_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include <optional>
#include "multisphere_datatypes.hpp"
#include "multisphere_io.hpp"
#include "multisphere_voxel_processing.hpp"

// --- Helper: Distance Squared ---
inline float get_dist_sq(const Eigen::Vector3f& a, const Eigen::Vector3f& b) {
    return (a - b).squaredNorm();
}

// -- Helper packed struct for cache-friendly sorting
struct PeakEntry {
    float val;      // Pre-fetched distance value
    int x, y, z;    // Coordinates
    float radius;   // Radius at this peak
    int index;
    // Original index (optional, if you needed to map back to original list)
};

/**
 * Detect local maxima in a 3D volume.
 * Returns an Nx4 matrix (Float) where columns are: x, y, z, radius
 */
Eigen::MatrixX4f peak_local_max_3d(const VoxelGrid<float>& field, const VoxelGrid<float>& original_distance, int min_distance = 2, int min_radius_vox = 1) {
    // CHANGE 1: Use Vector4f (floats)
    std::vector<Eigen::Vector4f> peaks;
    
    int nx = (int)field.nx(); 
    int ny = (int)field.ny();
    int nz = (int)field.nz();

    #pragma omp parallel
    {
        std::vector<Eigen::Vector4f> local_peaks;
        
        #pragma omp for collapse(2)
        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < ny; ++y) {
                for (int z = 0; z < nz; ++z) {
                    float val = field(x, y, z);
                    float radius = original_distance(x, y, z);
                    
                    if (radius < min_radius_vox) continue;
                    if (val <= 1.0f) continue;

                    bool is_max = true;
                    // Check local neighborhood
                    for (int dx = -min_distance; is_max && dx <= min_distance; ++dx) {
                        for (int dy = -min_distance; is_max && dy <= min_distance; ++dy) {
                            for (int dz = -min_distance; is_max && dz <= min_distance; ++dz) {
                                int nx_idx = x + dx, ny_idx = y + dy, nz_idx = z + dz;
                                
                                if (nx_idx >= 0 && nx_idx < nx && ny_idx >= 0 && ny_idx < ny && nz_idx >= 0 && nz_idx < nz) {
                                    if (field(nx_idx, ny_idx, nz_idx) > val) is_max = false;
                                }
                            }
                        }
                    }
                    
                    if (is_max) {
                        // CHANGE 2: Cast to float (lighter memory usage)
                        local_peaks.push_back(Eigen::Vector4f((float)x, (float)y, (float)z, radius));
                    }
                }
            }
        }
        #pragma omp critical
        peaks.insert(peaks.end(), local_peaks.begin(), local_peaks.end());
    }

    // CHANGE 3: Sort based on radius (descending), then x,y,z (ascending)
    std::sort(peaks.begin(), peaks.end(), [](const Eigen::Vector4f& a, const Eigen::Vector4f& b) {
        float val_a = a[3]; 
        float val_b = b[3];
        // Epsilon for float is usually larger than double, but 1e-5 is still fine for this scale
        const float epsilon = 1e-5f;
        if (std::abs(val_a - val_b) > epsilon) {
            return val_a > val_b;
        }
        if (a[0] != b[0]) return a[0] < b[0]; 
        if (a[1] != b[1]) return a[1] < b[1]; 
        return a[2] < b[2];                   
    });

    // CHANGE 4: Return MatrixX4f
    Eigen::MatrixX4f result(peaks.size(), 4);
    for (size_t i = 0; i < peaks.size(); ++i) {
        result.row(i) = peaks[i];
    }
    
    return result;
}


/**
 * Filters peaks based on minimum center distance and applies sub-voxel shifts.
 * Returns filtered peaks as an Nx4 matrix (Float).
 */

/**
 * Filters peaks based on minimum center distance against a history of spheres 
 * AND self-collisions, then applies sub-voxel shifts.
 * * Returns filtered peaks as an Nx4 matrix (x, y, z, radius).
 */
Eigen::MatrixX4f filter_and_shift_peaks(
    const Eigen::MatrixX4f& peaks,
    const Eigen::MatrixX4f& sphere_table, // Added: The History table
    VoxelGrid<float>& distance_field,
    int min_center_distance_vox) 
{
    // If no candidates, return immediately
    if (peaks.rows() == 0) return peaks;

    // Squared threshold for unified center-to-center distance
    const float min_dist_sq = (float)(min_center_distance_vox * min_center_distance_vox);

    // Storage for the final refined spheres
    std::vector<Eigen::Vector4f> accepted_spheres;
    accepted_spheres.reserve(peaks.rows());

    // --- The Scan ---
    for (int i = 0; i < peaks.rows(); ++i) {
        float px = peaks(i, 0);
        float py = peaks(i, 1);
        float pz = peaks(i, 2);
        
        bool collision = false;

        // 1. Check against ALL existing spheres (The History)
        // We do this first to discard peaks colliding with previously known structures
        
        for (int j = 0; j < sphere_table.rows(); ++j) {
            float hx = (float)sphere_table(j, 0);
            float hy = (float)sphere_table(j, 1);
            float hz = (float)sphere_table(j, 2);

            float dx = px - hx;
            float dy = py - hy;
            float dz = pz - hz;

            if ((dx*dx + dy*dy + dz*dz) < min_dist_sq) {
                collision = true;
                break;
            }
        }

        // If it collided with history, skip it; otherwise check self-collision
        if (!collision) {
            // 2. Check against ALREADY ACCEPTED peaks in this batch (Self-NMS)
            for (const auto& sphere : accepted_spheres) {
                float dx = px - sphere[0];
                float dy = py - sphere[1];
                float dz = pz - sphere[2];

                if ((dx*dx + dy*dy + dz*dz) < min_dist_sq) {
                    collision = true;
                    break;
                }
            }
        }

        // If valid (no collision with history OR new neighbors)
        if (!collision) {
            
            // Get integer coordinates for the voxel lookup
            Eigen::Vector3i coord_int = peaks.row(i).head<3>().cast<int>();

            // Perform the sub-voxel shift logic
            Eigen::Vector3f direction = shift_voxel_center(distance_field, coord_int, 0.25f);
            Eigen::Vector3f shift_vec = 0.5f * direction;
            Eigen::Vector3f true_center = coord_int.cast<float>() + shift_vec;

            // Calculate radius extending to voxel corners
            float true_radius = peaks(i, 3) + shift_vec.norm() ; 

            accepted_spheres.push_back(Eigen::Vector4f(
                true_center.x(), 
                true_center.y(), 
                true_center.z(), 
                true_radius
            ));
        }
    }

    // --- Pack Result ---
    Eigen::MatrixX4f result(accepted_spheres.size(), 4);
    for (size_t k = 0; k < accepted_spheres.size(); ++k) {
        result.row(k) = accepted_spheres[k];
    }

    return result;
}

// --- Append Sphere Table (In-Place) ---
// If max_spheres is set, respects that limit.
// If peaks.rows()+table.rows() < max_spheres, all peaks are added simultaneously.
void append_sphere_table(
    Eigen::MatrixX4f& table,                // Modified in place
    const Eigen::MatrixX4f& peaks,
    std::optional<int> max_spheres = std::nullopt) 
{
    if (peaks.rows() == 0) return;

    if (!max_spheres.has_value() || 
       (table.rows() + peaks.rows() < static_cast<size_t>(*max_spheres))) {
        table.conservativeResize(table.rows() + peaks.rows(), 4);
        for (int i = 0; i < peaks.rows(); ++i) {
            table.row(table.rows() - peaks.rows() + i) = peaks.row(i);
        }
        return;
    }
    else {
        // We can only add up to the max_spheres limit
        int space_left = *max_spheres - table.rows();
        if (space_left <= 0) return;

        table.conservativeResize(table.rows() + space_left, 4);
        for (int i = 0; i < space_left; ++i) {
            table.row(table.rows() - space_left + i) = peaks.row(i);
        }
    }
}

// --- Residual Distance Field ---
VoxelGrid<float> residual_distance_field(
    const VoxelGrid<float>& original_distance,
    const VoxelGrid<float>& spheres_distance
   ) 
{
    if (original_distance.nx() != spheres_distance.nx() || 
        original_distance.ny() != spheres_distance.ny() || 
        original_distance.nz() != spheres_distance.nz()) {
        throw std::invalid_argument("Shapes must match.");
    }

    VoxelGrid<float> residual(original_distance.nx(), original_distance.ny(), original_distance.nz(), 
                              original_distance.voxel_size, original_distance.origin);

    #pragma omp parallel for
    for (size_t i = 0; i < residual.data.size(); ++i) {
        // R = D_orig - D_spheres
        float val = original_distance.data[i] - spheres_distance.data[i];
        residual.data[i] = (val > 0.0f) ? val : 0.0f;
    }
    return residual;
}

#endif