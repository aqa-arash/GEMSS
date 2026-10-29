/**
 * @file GEMSS_edt.hpp
 * @brief Exact Euclidean distance transform of a binary voxel grid.
 *
 * Separable algorithm of Felzenszwalb & Huttenlocher, "Distance Transforms of Sampled Functions",
 * Theory of Computing 8 (2012) 415-428: a 1D distance pass along x on the binary mask, then the lower
 * envelope of parabolas along y and along z. Every intermediate value is a squared distance between voxel
 * centres, i.e. an integer (exact in float below 2^24), so the transform is exact; the result is its root.
 *
 * Conventions: index = (x * ny + y) * nz + z (z contiguous). Voxels > 0 are foreground, the rest background.
 * Everything outside the grid counts as background, one voxel beyond each face. A grid with nz == 1 is a
 * 2D slice: it has no border in z.
 *
 * @author Arash Moradian
 * @date 2026-10-29
 */

#ifndef GEMSS_EDT_HPP
#define GEMSS_EDT_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#ifdef HAVE_OPENMP
    #include <omp.h>
#endif

namespace GEMSS {
namespace edt_detail {

/// Per-thread scratch for one line of the lower envelope.
struct EnvelopeScratch {
    std::vector<double> f;  ///< samples of the line
    std::vector<int> v;     ///< apexes of the parabolas in the envelope
    std::vector<double> z;  ///< boundaries between envelope parabolas
    explicit EnvelopeScratch(int n) : f(n), v(n), z(n + 1) {}
};

/**
 * @brief Lower envelope along one line, in place: d[q] = min_p ((q - p)^2 + f[p]).
 * @param line First sample; samples are `stride` apart.
 * @param n Number of samples.
 * @param border Also take the background sites p = -1 and p = n into account.
 * @param take_root Write sqrt(d) instead of d.
 */
inline void envelope_1d(float* line, int n, long stride, bool border, bool take_root, EnvelopeScratch& s) {
    bool background_only = true;
    for (int q = 0; q < n; ++q) {
        s.f[q] = line[q * stride];
        background_only &= (s.f[q] == 0.0);
    }
    if (background_only) return;  // all zero, and zero is its own root

    const double inf = std::numeric_limits<double>::infinity();
    int k = 0;
    s.v[0] = 0;
    s.z[0] = -inf;
    s.z[1] = inf;
    for (int q = 1; q < n; ++q) {
        // intersection of the parabolas rooted at q and at the last envelope apex
        double x = ((s.f[q] + double(q) * q) - (s.f[s.v[k]] + double(s.v[k]) * s.v[k])) / (2.0 * (q - s.v[k]));
        while (x <= s.z[k]) {
            --k;
            x = ((s.f[q] + double(q) * q) - (s.f[s.v[k]] + double(s.v[k]) * s.v[k])) / (2.0 * (q - s.v[k]));
        }
        ++k;
        s.v[k] = q;
        s.z[k] = x;
        s.z[k + 1] = inf;
    }

    k = 0;
    for (int q = 0; q < n; ++q) {
        while (s.z[k + 1] < q) ++k;
        double d = double(q - s.v[k]) * (q - s.v[k]) + s.f[s.v[k]];
        if (border) d = std::min(d, std::min(double(q + 1) * (q + 1), double(n - q) * (n - q)));
        const float df = static_cast<float>(d);
        line[q * stride] = take_root ? std::sqrt(df) : df;
    }
}

} // namespace edt_detail

/**
 * @brief Exact Euclidean distance transform (in voxels) of a binary grid.
 * @param mask Input grid, nx * ny * nz values; > 0 is foreground.
 * @param out Output grid of the same size: distance of each voxel centre to the nearest background voxel centre.
 */
template <typename T>
inline void euclidean_distance_transform(const T* mask, int nx, int ny, int nz, float* out) {
    const long sx = static_cast<long>(ny) * nz;
    const bool border_z = nz > 1;

    #pragma omp parallel
    {
        // Pass 1, along x: distance to the nearest background voxel on the same x-line, squared.
        // Each thread sweeps whole xz-planes, so every inner loop runs over contiguous z.
        std::vector<float> run(nz);
        #pragma omp for schedule(static)
        for (int y = 0; y < ny; ++y) {
            std::fill(run.begin(), run.end(), 0.0f);  // background just before x = 0
            for (int x = 0; x < nx; ++x) {
                const T* m = mask + x * sx + static_cast<long>(y) * nz;
                float* o = out + x * sx + static_cast<long>(y) * nz;
                for (int z = 0; z < nz; ++z) {
                    run[z] = (m[z] > T(0)) ? run[z] + 1.0f : 0.0f;
                    o[z] = run[z];
                }
            }
            std::fill(run.begin(), run.end(), 0.0f);  // background just after x = nx - 1
            for (int x = nx - 1; x >= 0; --x) {
                float* o = out + x * sx + static_cast<long>(y) * nz;
                for (int z = 0; z < nz; ++z) {
                    const float d = std::min(o[z], run[z] + 1.0f);
                    run[z] = d;
                    o[z] = d * d;
                }
            }
        }

        // Pass 2, along y (stride nz), and pass 3, along z (contiguous, takes the root).
        edt_detail::EnvelopeScratch scratch(std::max(ny, nz));
        #pragma omp for collapse(2) schedule(dynamic, 16)
        for (int x = 0; x < nx; ++x)
            for (int z = 0; z < nz; ++z)
                edt_detail::envelope_1d(out + x * sx + z, ny, nz, true, false, scratch);

        #pragma omp for collapse(2) schedule(dynamic, 16)
        for (int x = 0; x < nx; ++x)
            for (int y = 0; y < ny; ++y)
                edt_detail::envelope_1d(out + x * sx + static_cast<long>(y) * nz, nz, 1, border_z, true, scratch);
    }
}

} // namespace GEMSS

#endif // GEMSS_EDT_HPP
