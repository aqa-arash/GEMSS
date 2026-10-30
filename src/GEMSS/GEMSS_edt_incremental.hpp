/**
 * @file GEMSS_edt_incremental.hpp
 * @brief Exact incremental Euclidean distance transform for a binary mask that only grows.
 *
 * The separable transform of GEMSS_edt.hpp runs three passes of 1D transforms:
 *   - pass 1 (A): squared distance along x, computed from the mask;
 *   - pass 2 (B): lower envelope of A along y;
 *   - pass 3 (C): lower envelope of B along z; the distance is sqrt(C).
 * Each line of a pass depends only on the same line of the previous pass. IncrementalEDT keeps A and B of the
 * last transform and, after voxels were added to the mask, recomputes only the lines whose input changed:
 *   - pass 1: the x-lines marked by the caller (lines with a voxel that changed from 0 to 1);
 *   - pass 2: the y-lines on which A changed;
 *   - pass 3: the z-lines on which B changed.
 * A recomputed line uses the same arithmetic as the full transform and a skipped line has an unchanged input,
 * so the result is bit-identical to a full transform. Every recomputed line is compared with its stored values,
 * so the update always reports where the distance changed.
 *
 * Define GEMSS_EDT_VERIFY to compare every update with a full transform (debugging only).
 *
 * @author Arash Moradian
 * @date 2026-10-29
 */

#ifndef GEMSS_EDT_INCREMENTAL_HPP
#define GEMSS_EDT_INCREMENTAL_HPP

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <vector>
#ifdef GEMSS_EDT_VERIFY
    #include <cstdio>
#endif
#include "GEMSS_datatypes.hpp"
#include "GEMSS_edt.hpp"

namespace GEMSS {

/**
 * @brief Incremental Euclidean distance transform of a growing binary mask.
 *
 * Mark the x-line (y, z) of every voxel that changes from 0 to 1 with mark_line(), then call update(). Any other
 * change of the mask requires invalidate(). After an update, changed_blocks() flags the blocks in which a
 * distance changed; the first update reports all_changed() instead.
 */
class IncrementalEDT {
public:
    /// If more than this fraction of the x-lines is dirty, pass 1 sweeps the whole grid plane by plane.
    double full_fraction = 0.5;

    /**
     * @brief Creates the transform for an nx x ny x nz grid. The first update is a full transform.
     * @param block Edge length (voxels) of the blocks reported by changed_blocks().
     */
    IncrementalEDT(int nx, int ny, int nz, int block = 8)
        : nx_(nx), ny_(ny), nz_(nz), stride_x_(static_cast<long>(ny) * nz),
          block_(block),
          block_grid_{(nx + block - 1) / block, (ny + block - 1) / block, (nz + block - 1) / block},
          changed_blocks_(static_cast<size_t>(block_grid_[0]) * block_grid_[1] * block_grid_[2], 0),
          a_(static_cast<size_t>(nx) * ny * nz),
          b_(static_cast<size_t>(nx) * ny * nz),
          dirty_x_(static_cast<size_t>(ny) * nz, 0),
          dirty_y_(static_cast<size_t>(nx) * nz, 0),
          dirty_z_(static_cast<size_t>(nx) * ny, 0) {}

    /// Marks x-line (y, z) as changed. Thread-safe.
    void mark_line(int y, int z) { set_flag(&dirty_x_[static_cast<size_t>(y) * nz_ + z]); }

    /// Forgets the stored passes, so that the next update is a full transform.
    void invalidate() { valid_ = false; }

    /// True if the last update was a complete transform (first update); changed_blocks() is then not filled.
    bool all_changed() const { return all_changed_; }

    /// Flags of the blocks in which the last update changed a distance, indexed (bx * nby + by) * nbz + bz.
    const std::vector<uint8_t>& changed_blocks() const { return changed_blocks_; }

    /// Edge length (voxels) of the blocks.
    int block_size() const { return block_; }

    /// Number of blocks along x, y and z.
    const std::array<int, 3>& block_grid() const { return block_grid_; }

    /**
     * @brief Updates the distance transform of `mask` in `dt`.
     * @param mask Binary mask (> 0 is foreground), same shape as given to the constructor.
     * @param dt Distance field; must hold the result of this object's previous update (any content on the first).
     * @throws std::invalid_argument if the mask shape does not match.
     */
    void update(const VoxelGrid<uint8_t>& mask, VoxelGrid<float>& dt) {
        if (static_cast<int>(mask.nx()) != nx_ || static_cast<int>(mask.ny()) != ny_ || static_cast<int>(mask.nz()) != nz_) {
            throw std::invalid_argument("IncrementalEDT: mask shape does not match.");
        }
        if (dt.shape != mask.shape) {
            dt = VoxelGrid<float>(mask.nx(), mask.ny(), mask.nz(), mask.voxel_size, mask.origin);
            valid_ = false;
        }
        std::fill(changed_blocks_.begin(), changed_blocks_.end(), 0);
        all_changed_ = false;

        const uint8_t* m = mask.data.data();
        float* c = dt.data.data();
        if (valid_) {
            update_lines(m, c);
        } else {
            // Stored values may be stale: recompute every line (writing only differing values is still complete).
            full_pass_x(m);
            update_y_lines(all_lines(dirty_y_.size()));
            update_z_lines(all_lines(dirty_z_.size()), c);
            all_changed_ = true;
            valid_ = true;
        }
        std::fill(dirty_x_.begin(), dirty_x_.end(), 0);

        #ifdef GEMSS_EDT_VERIFY
        verify(mask, c);
        #endif
    }

private:
    int nx_, ny_, nz_;
    long stride_x_;                       ///< index step of x (z is contiguous, y has step nz)
    bool valid_ = false;                  ///< a_, b_ and the caller's distance field hold the last transform
    int block_;                           ///< block edge length (voxels)
    std::array<int, 3> block_grid_;       ///< number of blocks along x, y, z
    bool all_changed_ = false;            ///< the last update was a complete transform
    std::vector<uint8_t> changed_blocks_; ///< blocks in which the last update changed a distance
    std::vector<float> a_;                ///< pass-1 output (squared distance along x)
    std::vector<float> b_;                ///< pass-2 output (squared distance in the xy-slice)
    std::vector<uint8_t> dirty_x_;        ///< x-lines to recompute, indexed y * nz + z
    std::vector<uint8_t> dirty_y_;        ///< y-lines to recompute, indexed x * nz + z
    std::vector<uint8_t> dirty_z_;        ///< z-lines to recompute, indexed x * ny + y

    /// Sets a flag to 1 from any thread.
    static void set_flag(uint8_t* flag) { __atomic_store_n(flag, static_cast<uint8_t>(1), __ATOMIC_RELAXED); }

    /// Indices of the set flags, in increasing order.
    static std::vector<int> set_indices(const std::vector<uint8_t>& flags) {
        std::vector<int> indices;
        for (size_t i = 0; i < flags.size(); ++i) {
            if (flags[i]) indices.push_back(static_cast<int>(i));
        }
        return indices;
    }

    /// Indices 0 .. count-1.
    static std::vector<int> all_lines(size_t count) {
        std::vector<int> indices(count);
        for (size_t i = 0; i < count; ++i) indices[i] = static_cast<int>(i);
        return indices;
    }

    /// Pass 1 on every x-line, sweeping the grid plane by plane: a_ = squared distance along x.
    void full_pass_x(const uint8_t* m) {
        const int nx = nx_, ny = ny_, nz = nz_;
        const long sx = stride_x_;
        #pragma omp parallel
        {
            std::vector<float> run(nz);
            #pragma omp for schedule(static)
            for (int y = 0; y < ny; ++y) {
                std::fill(run.begin(), run.end(), 0.0f);  // background just before x = 0
                for (int x = 0; x < nx; ++x) {
                    const uint8_t* __restrict mm = m + x * sx + static_cast<long>(y) * nz;
                    float* __restrict out = a_.data() + x * sx + static_cast<long>(y) * nz;
                    float* __restrict r = run.data();
                    for (int z = 0; z < nz; ++z) {
                        const int foreground = mm[z] > 0;  // branch-free so the loop vectorizes
                        r[z] = (r[z] + 1.0f) * static_cast<float>(foreground);
                        out[z] = r[z];
                    }
                }
                std::fill(run.begin(), run.end(), 0.0f);  // background just after x = nx - 1
                for (int x = nx - 1; x >= 0; --x) {
                    float* out = a_.data() + x * sx + static_cast<long>(y) * nz;
                    for (int z = 0; z < nz; ++z) {
                        const float d = std::min(out[z], run[z] + 1.0f);
                        run[z] = d;
                        out[z] = d * d;
                    }
                }
            }
        }
    }

    // ---- Line-by-line update ----

    /**
     * @brief Recomputes the lines whose input changed, pass by pass.
     *
     * Passes 2 and 3 always run on exactly the dirty lines: when most lines are dirty, the list is in grid order
     * and costs about as much as a full pass, and the comparison still reports what changed. Pass 1 sweeps the
     * whole grid when many x-lines are dirty, since single x-lines are strided in memory; every y-line then
     * counts as dirty.
     */
    void update_lines(const uint8_t* m, float* c) {
        const std::vector<int> x_lines = set_indices(dirty_x_);
        if (x_lines.empty()) return;  // mask unchanged
        if (static_cast<double>(x_lines.size()) > full_fraction * dirty_x_.size()) {
            full_pass_x(m);
            std::fill(dirty_y_.begin(), dirty_y_.end(), 1);
        } else {
            update_x_lines(m, x_lines);
        }

        const std::vector<int> y_lines = set_indices(dirty_y_);
        if (y_lines.empty()) return;
        update_y_lines(y_lines);

        const std::vector<int> z_lines = set_indices(dirty_z_);
        if (z_lines.empty()) return;
        update_z_lines(z_lines, c);
    }

    /// Pass 1 on the given x-lines (index y * nz + z); flags the y-lines on which a_ changed.
    void update_x_lines(const uint8_t* m, const std::vector<int>& lines) {
        const int nx = nx_, nz = nz_;
        const long sx = stride_x_;
        std::fill(dirty_y_.begin(), dirty_y_.end(), 0);
        #pragma omp parallel
        {
            std::vector<float> line(nx);
            #pragma omp for schedule(dynamic, 64)
            for (size_t l = 0; l < lines.size(); ++l) {
                const int y = lines[l] / nz, z = lines[l] % nz;
                const long base = static_cast<long>(y) * nz + z;
                float run = 0.0f;
                for (int x = 0; x < nx; ++x) {
                    const int foreground = m[x * sx + base] > 0;
                    run = (run + 1.0f) * static_cast<float>(foreground);
                    line[x] = run;
                }
                run = 0.0f;
                for (int x = nx - 1; x >= 0; --x) {
                    const float d = std::min(line[x], run + 1.0f);
                    run = d;
                    line[x] = d * d;
                }
                for (int x = 0; x < nx; ++x) {
                    float& a = a_[x * sx + base];
                    if (line[x] != a) {
                        a = line[x];
                        set_flag(&dirty_y_[static_cast<long>(x) * nz + z]);
                    }
                }
            }
        }
    }

    /// Pass 2 on the given y-lines (index x * nz + z); flags the z-lines on which b_ changed.
    void update_y_lines(const std::vector<int>& lines) {
        const int ny = ny_, nz = nz_;
        const long sx = stride_x_;
        std::fill(dirty_z_.begin(), dirty_z_.end(), 0);
        #pragma omp parallel
        {
            edt_detail::EnvelopeScratch scratch(std::max(ny, nz));
            std::vector<float> line(ny);
            #pragma omp for schedule(dynamic, 16)
            for (size_t l = 0; l < lines.size(); ++l) {
                const int x = lines[l] / nz, z = lines[l] % nz;
                const long base = x * sx + z;
                for (int y = 0; y < ny; ++y) line[y] = a_[base + static_cast<long>(y) * nz];
                edt_detail::envelope_1d(line.data(), ny, 1, true, false, scratch);
                for (int y = 0; y < ny; ++y) {
                    float& b = b_[base + static_cast<long>(y) * nz];
                    if (line[y] != b) {
                        b = line[y];
                        set_flag(&dirty_z_[static_cast<long>(x) * ny + y]);
                    }
                }
            }
        }
    }

    /// Pass 3 on the given z-lines (index x * ny + y), writing the distance; flags the blocks in which it changed.
    void update_z_lines(const std::vector<int>& lines, float* c) {
        const int ny = ny_, nz = nz_;
        const long sx = stride_x_;
        #pragma omp parallel
        {
            edt_detail::EnvelopeScratch scratch(std::max(ny, nz));
            std::vector<float> line(nz);
            #pragma omp for schedule(dynamic, 16)
            for (size_t l = 0; l < lines.size(); ++l) {
                const int x = lines[l] / ny, y = lines[l] % ny;
                const long base = x * sx + static_cast<long>(y) * nz;
                std::memcpy(line.data(), b_.data() + base, nz * sizeof(float));
                edt_detail::envelope_1d(line.data(), nz, 1, nz > 1, true, scratch);
                const long block_row = (static_cast<long>(x / block_) * block_grid_[1] + y / block_) * block_grid_[2];
                for (int z = 0; z < nz; ++z) {
                    if (line[z] != c[base + z]) {
                        c[base + z] = line[z];
                        set_flag(&changed_blocks_[block_row + z / block_]);
                    }
                }
            }
        }
    }

    #ifdef GEMSS_EDT_VERIFY
    /// Compares `c` with a full transform of `mask` and reports differing voxels on stderr.
    void verify(const VoxelGrid<uint8_t>& mask, const float* c) const {
        std::vector<float> reference(mask.data.size());
        euclidean_distance_transform(mask.data.data(), nx_, ny_, nz_, reference.data());
        size_t differing = 0;
        for (size_t i = 0; i < reference.size(); ++i) {
            differing += std::memcmp(&reference[i], &c[i], sizeof(float)) != 0;
        }
        if (differing > 0) {
            std::fprintf(stderr, "[GEMSS_EDT_VERIFY] incremental EDT differs from full transform in %zu voxels\n", differing);
        }
    }
    #endif
};

} // namespace GEMSS

#endif // GEMSS_EDT_INCREMENTAL_HPP
