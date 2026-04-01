#ifndef GEMSS_PHYSICS_COMPUTATION_HPP
#define GEMSS_PHYSICS_COMPUTATION_HPP

/**
 * @file multisphere_physics_computation.hpp
 * @brief Physical property calculations for GEMSS.
 *
 * Computes volume, center of mass, and inertia tensor from a voxel grid mask
 * using a highly-optimized, single-pass Parallel Axis Theorem reduction.
 *
 * @author Arash Moradian
 * @date 2026-03-12
 */

#include <cmath>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#ifdef HAVE_OPENMP
    #include <omp.h>
#endif

#include "GEMSS_datatypes.hpp"

namespace GEMSS {

/**
 * @brief Computes volume, Center of Mass (CoM), and inertia tensor from the final voxel mask in a single pass.
 * Updates the physical properties directly inside the provided SpherePack.
 * @param pack The SpherePack object to update.
 * @param voxelGrid The binary voxel grid.
 */
inline void compute_multisphere_physics(SpherePack& pack, const VoxelGrid<uint8_t>& voxelGrid) {
    long long N_vox = 0;
    
    // Moments accumulated in LOCAL physical space (origin subtracted)
    double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
    double sum_xx = 0.0, sum_yy = 0.0, sum_zz = 0.0;
    double sum_xy = 0.0, sum_xz = 0.0, sum_yz = 0.0;

    int nx = voxelGrid.nx();
    int ny = voxelGrid.ny();
    int nz = voxelGrid.nz();
    double vs = voxelGrid.voxel_size;

    #pragma omp parallel for collapse(3) reduction(+:N_vox, sum_x, sum_y, sum_z, sum_xx, sum_yy, sum_zz, sum_xy, sum_xz, sum_yz)
    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int z = 0; z < nz; ++z) {
                if (voxelGrid(x, y, z) > 0) {
                    N_vox++;
                    
                    // Local physical coordinates
                    double cx = x * vs;
                    double cy = y * vs;
                    double cz = z * vs;

                    sum_x += cx;
                    sum_y += cy;
                    sum_z += cz;
                    
                    sum_xx += cx * cx;
                    sum_yy += cy * cy;
                    sum_zz += cz * cz;
                    
                    sum_xy += cx * cy;
                    sum_xz += cx * cz;
                    sum_yz += cy * cz;
                }
            }
        }
    }

    if (N_vox == 0) {
        pack.volume = 0.0;
        pack.center_of_mass.setZero();
        pack.inertia_tensor.setZero();
        pack.principal_moments.setZero();
        pack.principal_axes.setIdentity();
        return;
    }

    double voxel_vol = vs * vs * vs;
    pack.volume = N_vox * voxel_vol;

    // 1. Calculate LOCAL Center of Mass
    double cx_local = sum_x / N_vox;
    double cy_local = sum_y / N_vox;
    double cz_local = sum_z / N_vox;

    // 2. Shift Local CoM to GLOBAL CoM
    pack.center_of_mass << cx_local + voxelGrid.origin.x(),
                           cy_local + voxelGrid.origin.y(),
                           cz_local + voxelGrid.origin.z();

    // 3. Parallel Axis Theorem (Calculated entirely in stable local space)
    double Ixx = (sum_yy * voxel_vol) + (sum_zz * voxel_vol) - pack.volume * (cy_local * cy_local + cz_local * cz_local);
    double Iyy = (sum_xx * voxel_vol) + (sum_zz * voxel_vol) - pack.volume * (cx_local * cx_local + cz_local * cz_local);
    double Izz = (sum_xx * voxel_vol) + (sum_yy * voxel_vol) - pack.volume * (cx_local * cx_local + cy_local * cy_local);

    double Ixy = -(sum_xy * voxel_vol) + (pack.volume * cx_local * cy_local);
    double Ixz = -(sum_xz * voxel_vol) + (pack.volume * cx_local * cz_local);
    double Iyz = -(sum_yz * voxel_vol) + (pack.volume * cy_local * cz_local);

    double self_inertia_total = pack.volume * (vs * vs) / 6.0;

    pack.inertia_tensor << Ixx + self_inertia_total, Ixy, Ixz,
                           Ixy, Iyy + self_inertia_total, Iyz,
                           Ixz, Iyz, Izz + self_inertia_total;

    // 4. Principal Moments and Axes with Analytical Solver (Highly optimized for 3x3 symmetric tensors)                       
    // 4.1. Pre-condition the Inertia Tensor (Eliminates noise causing arbitrary degenerate rotations)
    // Use a threshold relative to the largest principal moment
    float max_I = pack.inertia_tensor.diagonal().maxCoeff();
    const float noise_threshold = max_I * 1e-5f; 

    pack.inertia_tensor = pack.inertia_tensor.unaryExpr([noise_threshold](float v) {
        return std::abs(v) < noise_threshold ? 0.0f : v;
    });

    // 4.2. Fast Analytical 3x3 Solver (Avoids iterative overhead)
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigensolver;
    eigensolver.computeDirect(pack.inertia_tensor);

    if (eigensolver.info() == Eigen::Success) {
        Eigen::Matrix3f evecs = eigensolver.eigenvectors();

        // 4.3. deterministic sign and handedness enforcement
        for (int i = 0; i < 3; ++i) {
            int max_idx;
            evecs.col(i).cwiseAbs().maxCoeff(&max_idx);
            // Flip if the dominant component is negative
            if (evecs(max_idx, i) < 0.0f) {
                evecs.col(i) *= -1.0f;
            }
        }

        // Force right-handed coordinate system
        if (evecs.col(0).cross(evecs.col(1)).dot(evecs.col(2)) < 0.0f) {
            evecs.col(2) *= -1.0f;
        }

        pack.principal_moments = eigensolver.eigenvalues();
        pack.principal_axes = evecs;
    } else {
        std::cerr << "[Warning] Analytical eigendecomposition failed." << std::endl;
        pack.principal_moments.setZero();
        pack.principal_axes.setIdentity();
    }

    // 5. Calculate Bounding Radius
    double max_r = 0.0;
    for (size_t i = 0; i < pack.num_spheres(); ++i) {
        Eigen::Vector3f center = pack.centers.row(i);
        float radius = pack.radii(i);

        // Ensure pack.center_of_mass is cast to Vector3f for type consistency
        double dist = (center - pack.center_of_mass.cast<float>()).norm() + radius + std::sqrt(3) * vs; // Add a diagonal voxel for safety
        if (dist > max_r) max_r = dist;
    }
    pack.bounding_radius = static_cast<float>(max_r);
}

} // namespace MSS
#endif // GEMSS_PHYSICS_COMPUTATION_HPP
