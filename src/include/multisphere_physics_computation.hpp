#ifndef MULTISPHERE_PHYSICS_COMPUTATION_HPP
#define MULTISPHERE_PHYSICS_COMPUTATION_HPP

/**
 * @file multisphere_physics_computation.hpp
 * @brief Physical property calculations for multisphere-cpp.
 *
 * Computes volume, center of mass, and inertia tensor from a voxel grid mask
 * using a highly-optimized, single-pass Parallel Axis Theorem reduction.
 *
 * @author Arash Moradian
 * @date 2026-03-12
 */

#include <cmath>
#include <iostream>
#include <Eigen/Dense>

#ifdef HAVE_OPENMP
    #include <omp.h>
#endif

#include "multisphere_datatypes.hpp"

namespace MSS {

/**
 * @brief Computes volume, Center of Mass (CoM), and inertia tensor from the final voxel mask in a single pass.
 * Updates the physical properties directly inside the provided SpherePack.
 * @param pack The SpherePack object to update.
 * @param final_mask The binary voxel grid representing the multisphere union.
 */
inline void compute_multisphere_physics(SpherePack& pack, const VoxelGrid<uint8_t>& final_mask) {
    long long N_vox = 0;
    
    // Moments accumulated in LOCAL physical space (origin subtracted)
    double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
    double sum_xx = 0.0, sum_yy = 0.0, sum_zz = 0.0;
    double sum_xy = 0.0, sum_xz = 0.0, sum_yz = 0.0;

    int nx = final_mask.nx();
    int ny = final_mask.ny();
    int nz = final_mask.nz();
    double vs = final_mask.voxel_size;

    #pragma omp parallel for collapse(3) reduction(+:N_vox, sum_x, sum_y, sum_z, sum_xx, sum_yy, sum_zz, sum_xy, sum_xz, sum_yz)
    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int z = 0; z < nz; ++z) {
                if (final_mask(x, y, z) > 0) {
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
    pack.center_of_mass << cx_local + final_mask.origin.x(),
                           cy_local + final_mask.origin.y(),
                           cz_local + final_mask.origin.z();

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

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigensolver(pack.inertia_tensor);
    if (eigensolver.info() == Eigen::Success) {
        pack.principal_moments = eigensolver.eigenvalues();
        pack.principal_axes = eigensolver.eigenvectors();
    } else {
        std::cerr << "[Warning] Eigendecomposition failed." << std::endl;
        pack.principal_moments.setZero();
        pack.principal_axes.setIdentity();
    }
}

} // namespace MSS

#endif // MULTISPHERE_PHYSICS_COMPUTATION_HPP