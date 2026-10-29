/**
 * @file GEMSS_edt_incremental.hpp
 * @brief Exact incremental Euclidean distance transform for a binary mask that only grows.
 *
 * The separable transform of GEMSS_edt.hpp computes
 *   A: pass 1, squared 1D distance along x, from the mask,
 *   B: pass 2, lower envelope of A along y,
 *   C: pass 3, lower envelope of B along z (the result is sqrt(C)).
 * Every line of a pass depends only on the same line of the previous pass. This class keeps A and B of the
 * last transform and, after voxels were added to the mask, recomputes only
 *   pass 1: the x-lines the caller marked (lines containing a voxel that changed from 0 to 1),
 *   pass 2: the y-lines on which A changed,
 *   pass 3: the z-lines on which B changed.
 * A recomputed line uses the same arithmetic as the full transform and a skipped line has an unchanged input,
 * so the result is bit-identical to a full transform. A pass with more than `full_fraction` dirty lines is
 * recomputed in full (also exact).
 *
 * After each update, all_changed() and changed_blocks() tell which blocks (block^3 voxels) contain a voxel whose
 * distance changed, so later stages can restrict their work to those regions.
 *
 * Contract: between two updates the mask may only change by setting voxels from 0 to non-zero, and the x-line
 * (y, z) of every such voxel must be passed to mark_line(); otherwise call invalidate().
 * Define GEMSS_EDT_VERIFY to compare every update with a full transform.
 *
 * @author Arash Moradian
 * @date 2026-10-29
 */

#ifndef GEMSS_EDT_INCREMENTAL_HPP
#define GEMSS_EDT_INCREMENTAL_HPP

#include <algorithm>
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

class IncrementalEDT {
public:
    double full_fraction = 0.5;  ///< a pass with a larger fraction of dirty lines is recomputed in full

    IncrementalEDT(int nx, int ny, int nz, int block = 8)
        : nx_(nx), ny_(ny), nz_(nz), sx_(static_cast<long>(ny) * nz),
          block_(block), nbx_((nx + block - 1) / block), nby_((ny + block - 1) / block), nbz_((nz + block - 1) / block),
          changed_blocks_(static_cast<size_t>(nbx_) * nby_ * nbz_, 0),
          a_(static_cast<size_t>(nx) * ny * nz), b_(static_cast<size_t>(nx) * ny * nz),
          dirty1_(static_cast<size_t>(ny) * nz, 0), dirty2_(static_cast<size_t>(nx) * nz, 0),
          dirty3_(static_cast<size_t>(nx) * ny, 0) {}

    /// Mark x-line (y, z) as changed. Thread-safe.
    void mark_line(int y, int z) { set_flag(&dirty1_[static_cast<size_t>(y) * nz_ + z]); }

    /// Forget the stored passes; the next update is a full transform.
    void invalidate() { valid_ = false; }

    /// After update(): the whole field may have changed (a full pass ran); changed_blocks() is then not filled.
    bool all_changed() const { return all_changed_; }
    /// After update(): flags of the blocks containing a voxel whose distance changed, (bx * nby + by) * nbz + bz.
    const std::vector<uint8_t>& changed_blocks() const { return changed_blocks_; }
    int block_size() const { return block_; }

    /// Distance transform of `mask` into `dt`. `dt` must hold the result of this object's previous update.
    void update(const VoxelGrid<uint8_t>& mask, VoxelGrid<float>& dt) {
        if (static_cast<int>(mask.nx()) != nx_ || static_cast<int>(mask.ny()) != ny_ || static_cast<int>(mask.nz()) != nz_)
            throw std::invalid_argument("IncrementalEDT: mask shape does not match.");
        if (dt.shape != mask.shape) { dt = VoxelGrid<float>(mask.nx(), mask.ny(), mask.nz(), mask.voxel_size, mask.origin); valid_ = false; }
        const uint8_t* m = mask.data.data();
        float* C = dt.data.data();
        std::fill(changed_blocks_.begin(), changed_blocks_.end(), 0);
        all_changed_ = false;
        if (valid_) incremental(m, C);
        else { pass1_full(m); pass2_full(); pass3_full(C); valid_ = true; }
        std::fill(dirty1_.begin(), dirty1_.end(), 0);
        #ifdef GEMSS_EDT_VERIFY
        verify(mask, C);
        #endif
    }

private:
    int nx_, ny_, nz_;
    long sx_;
    bool valid_ = false;
    int block_, nbx_, nby_, nbz_;
    bool all_changed_ = false;
    std::vector<uint8_t> changed_blocks_;           ///< blocks with a changed distance (last update)
    std::vector<float> a_, b_;                       ///< pass-1 and pass-2 outputs of the last update
    std::vector<uint8_t> dirty1_, dirty2_, dirty3_;  ///< dirty x-lines (y,z), y-lines (x,z), z-lines (x,y)

    static void set_flag(uint8_t* flag) { __atomic_store_n(flag, static_cast<uint8_t>(1), __ATOMIC_RELAXED); }

    static std::vector<int> collect(const std::vector<uint8_t>& flags) {
        std::vector<int> list;
        for (size_t i = 0; i < flags.size(); ++i)
            if (flags[i]) list.push_back(static_cast<int>(i));
        return list;
    }

    // ---- full passes: same arithmetic as euclidean_distance_transform ----

    void pass1_full(const uint8_t* m) {
        const int nx = nx_, ny = ny_, nz = nz_;
        const long sx = sx_;
        #pragma omp parallel
        {
            std::vector<float> run(nz);
            #pragma omp for schedule(static)
            for (int y = 0; y < ny; ++y) {
                std::fill(run.begin(), run.end(), 0.0f);
                for (int x = 0; x < nx; ++x) {
                    const uint8_t* __restrict mm = m + x * sx + static_cast<long>(y) * nz;
                    float* __restrict o = a_.data() + x * sx + static_cast<long>(y) * nz;
                    float* __restrict r = run.data();
                    for (int z = 0; z < nz; ++z) {
                        const int foreground = mm[z] > 0;
                        r[z] = (r[z] + 1.0f) * static_cast<float>(foreground);
                        o[z] = r[z];
                    }
                }
                std::fill(run.begin(), run.end(), 0.0f);
                for (int x = nx - 1; x >= 0; --x) {
                    float* o = a_.data() + x * sx + static_cast<long>(y) * nz;
                    for (int z = 0; z < nz; ++z) {
                        const float d = std::min(o[z], run[z] + 1.0f);
                        run[z] = d;
                        o[z] = d * d;
                    }
                }
            }
        }
    }

    void pass2_full() {
        const int nx = nx_, ny = ny_, nz = nz_;
        const long sx = sx_;
        std::memcpy(b_.data(), a_.data(), a_.size() * sizeof(float));
        #pragma omp parallel
        {
            edt_detail::EnvelopeScratch s(std::max(ny, nz));
            #pragma omp for collapse(2) schedule(dynamic, 16)
            for (int x = 0; x < nx; ++x)
                for (int z = 0; z < nz; ++z)
                    edt_detail::envelope_1d(b_.data() + x * sx + z, ny, nz, true, false, s);
        }
    }

    void pass3_full(float* C) {
        all_changed_ = true;
        const int nx = nx_, ny = ny_, nz = nz_;
        const long sx = sx_;
        std::memcpy(C, b_.data(), b_.size() * sizeof(float));
        #pragma omp parallel
        {
            edt_detail::EnvelopeScratch s(std::max(ny, nz));
            #pragma omp for collapse(2) schedule(dynamic, 16)
            for (int x = 0; x < nx; ++x)
                for (int y = 0; y < ny; ++y)
                    edt_detail::envelope_1d(C + x * sx + static_cast<long>(y) * nz, nz, 1, nz > 1, true, s);
        }
    }

    // ---- incremental update ----

    void incremental(const uint8_t* m, float* C) {
        const int nx = nx_, ny = ny_, nz = nz_;
        const long sx = sx_;

        const std::vector<int> L1 = collect(dirty1_);
        if (L1.empty()) return;  // mask unchanged: C is still valid
        if (static_cast<double>(L1.size()) > full_fraction * dirty1_.size()) {
            pass1_full(m); pass2_full(); pass3_full(C);
            return;
        }

        // Pass 1 on the marked x-lines; mark the y-lines (x,z) where A changed.
        std::fill(dirty2_.begin(), dirty2_.end(), 0);
        #pragma omp parallel
        {
            std::vector<float> line(nx);
            #pragma omp for schedule(dynamic, 64)
            for (size_t l = 0; l < L1.size(); ++l) {
                const int y = L1[l] / nz, z = L1[l] % nz;
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
                    if (line[x] != a) { a = line[x]; set_flag(&dirty2_[static_cast<long>(x) * nz + z]); }
                }
            }
        }
        const std::vector<int> L2 = collect(dirty2_);
        if (L2.empty()) return;
        if (static_cast<double>(L2.size()) > full_fraction * dirty2_.size()) {
            pass2_full(); pass3_full(C);
            return;
        }

        // Pass 2 on the dirty y-lines; mark the z-lines (x,y) where B changed.
        std::fill(dirty3_.begin(), dirty3_.end(), 0);
        #pragma omp parallel
        {
            edt_detail::EnvelopeScratch s(std::max(ny, nz));
            std::vector<float> line(ny);
            #pragma omp for schedule(dynamic, 16)
            for (size_t l = 0; l < L2.size(); ++l) {
                const int x = L2[l] / nz, z = L2[l] % nz;
                const long base = x * sx + z;
                for (int y = 0; y < ny; ++y) line[y] = a_[base + static_cast<long>(y) * nz];
                edt_detail::envelope_1d(line.data(), ny, 1, true, false, s);
                for (int y = 0; y < ny; ++y) {
                    float& b = b_[base + static_cast<long>(y) * nz];
                    if (line[y] != b) { b = line[y]; set_flag(&dirty3_[static_cast<long>(x) * ny + y]); }
                }
            }
        }
        const std::vector<int> L3 = collect(dirty3_);
        if (L3.empty()) return;
        if (static_cast<double>(L3.size()) > full_fraction * dirty3_.size()) {
            pass3_full(C);
            return;
        }

        // Pass 3 on the dirty z-lines (contiguous), with the root; flag the blocks where C changed.
        #pragma omp parallel
        {
            edt_detail::EnvelopeScratch s(std::max(ny, nz));
            std::vector<float> line(nz);
            #pragma omp for schedule(dynamic, 16)
            for (size_t l = 0; l < L3.size(); ++l) {
                const int x = L3[l] / ny, y = L3[l] % ny;
                float* c = C + x * sx + static_cast<long>(y) * nz;
                std::memcpy(line.data(), b_.data() + x * sx + static_cast<long>(y) * nz, nz * sizeof(float));
                edt_detail::envelope_1d(line.data(), nz, 1, nz > 1, true, s);
                const long row_block = (static_cast<long>(x / block_) * nby_ + y / block_) * nbz_;
                for (int z = 0; z < nz; ++z) {
                    if (line[z] != c[z]) { c[z] = line[z]; set_flag(&changed_blocks_[row_block + z / block_]); }
                }
            }
        }
    }

    #ifdef GEMSS_EDT_VERIFY
    void verify(const VoxelGrid<uint8_t>& mask, const float* C) const {
        std::vector<float> ref(mask.data.size());
        euclidean_distance_transform(mask.data.data(), nx_, ny_, nz_, ref.data());
        size_t bad = 0;
        for (size_t i = 0; i < ref.size(); ++i) bad += std::memcmp(&ref[i], &C[i], sizeof(float)) != 0;
        if (bad) std::fprintf(stderr, "[GEMSS_EDT_VERIFY] incremental EDT differs from full transform in %zu voxels\n", bad);
    }
    #endif
};

} // namespace GEMSS

#endif // GEMSS_EDT_INCREMENTAL_HPP
