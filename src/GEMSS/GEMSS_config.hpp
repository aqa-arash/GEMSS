#ifndef GEMSS_CONFIG_HPP
#define GEMSS_CONFIG_HPP


#include <Eigen/Core>
#include <istream>
#include <limits>
#include <ostream>
#include <string>

namespace GEMSS {


struct MultisphereConfig {
    // --- Required for mesh reconstruction ---
    int div = 100;               // grid resolution
    int padding = 2;             // grid padding
    float minimum_radius_real = 0.0f; // min sphere radius (real units)
    bool confine_mesh = false;   // confine to mesh boundary

    // --- Peak detection and filtering ---
    float min_center_distance_rel = 1.4f; // min center distance (relative) computed as k <= sqrt(min_radius)
    int search_window = 2;       // neighbor distance for peaks
    float radius_offset_vox = 0.5f; // radius correction
    int min_radius_vox = 2;      // min sphere radius (voxels)


    // --- Convergence & limits ---
    float precision_target = 1.0f; // target precision
    int max_spheres = 0;         // max spheres (0 = unlimited)
    int persistence = 2;         // solver persistence

    // --- Utilities & prior state ---
    int compute_physics = 2;     // compute physics properties
    bool is_2d = false;   // treat reconstruction as planar (single z-slice + 2D physics)
    float density = 1.0f;        // material density for mass properties 
    bool prune_isolated_spheres = false; // remove isolated spheres
    bool show_progress = true;   // print progress
    Eigen::MatrixX4f initial_sphere_table = Eigen::MatrixX4f(0, 4); // initial spheres
    
};


// --- Text I/O (restart files) ---------------------------------------------------------------
// Every setting as "name value" pairs in declaration order, framed by "MultisphereConfig ... end", so a reader
// can check each name and a missing or renamed setting fails the read (failbit) instead of shifting values.
// initial_sphere_table (prior state of a reconstruction) is not written; a read leaves it as it was.

inline std::ostream& operator<<(std::ostream& os, const MultisphereConfig& c)
{
    const std::streamsize precision = os.precision(std::numeric_limits<float>::max_digits10);
    os << "MultisphereConfig"
       << " div " << c.div
       << " volumetric_div " << c.volumetric_div
       << " padding " << c.padding
       << " minimum_radius_real " << c.minimum_radius_real
       << " confine_mesh " << c.confine_mesh
       << " min_center_distance_rel " << c.min_center_distance_rel
       << " search_window " << c.search_window
       << " radius_offset_vox " << c.radius_offset_vox
       << " min_radius_vox " << c.min_radius_vox
       << " medial_flux_threshold " << c.medial_flux_threshold
       << " precision_target " << c.precision_target
       << " max_spheres " << c.max_spheres
       << " persistence " << c.persistence
       << " compute_physics " << c.compute_physics
       << " is_2d " << c.is_2d
       << " density " << c.density
       << " prune_isolated_spheres " << c.prune_isolated_spheres
       << " show_progress " << c.show_progress
       << " conserve_mass " << c.conserve_mass
       << " parent_representation_offset " << c.parent_representation_offset
       << " fracture_plane_offset " << c.fracture_plane_offset
       << " end";
    os.precision(precision);
    return os;
}

namespace detail {
// Reads "name value" and checks the name; on a mismatch the stream fails and the value is left unchanged
template <typename T>
inline void readNamed(std::istream& is, const char* name, T& value)
{
    std::string key;
    if (!(is >> key) || key != name) { is.setstate(std::ios::failbit); return; }
    is >> value;
}
} // namespace detail

inline std::istream& operator>>(std::istream& is, MultisphereConfig& c)
{
    std::string word;
    if (!(is >> word) || word != "MultisphereConfig") { is.setstate(std::ios::failbit); return is; }
    detail::readNamed(is, "div", c.div);
    detail::readNamed(is, "volumetric_div", c.volumetric_div);
    detail::readNamed(is, "padding", c.padding);
    detail::readNamed(is, "minimum_radius_real", c.minimum_radius_real);
    detail::readNamed(is, "confine_mesh", c.confine_mesh);
    detail::readNamed(is, "min_center_distance_rel", c.min_center_distance_rel);
    detail::readNamed(is, "search_window", c.search_window);
    detail::readNamed(is, "radius_offset_vox", c.radius_offset_vox);
    detail::readNamed(is, "min_radius_vox", c.min_radius_vox);
    detail::readNamed(is, "medial_flux_threshold", c.medial_flux_threshold);
    detail::readNamed(is, "precision_target", c.precision_target);
    detail::readNamed(is, "max_spheres", c.max_spheres);
    detail::readNamed(is, "persistence", c.persistence);
    detail::readNamed(is, "compute_physics", c.compute_physics);
    detail::readNamed(is, "is_2d", c.is_2d);
    detail::readNamed(is, "density", c.density);
    detail::readNamed(is, "prune_isolated_spheres", c.prune_isolated_spheres);
    detail::readNamed(is, "show_progress", c.show_progress);
    detail::readNamed(is, "conserve_mass", c.conserve_mass);
    detail::readNamed(is, "parent_representation_offset", c.parent_representation_offset);
    detail::readNamed(is, "fracture_plane_offset", c.fracture_plane_offset);
    if (!(is >> word) || word != "end") is.setstate(std::ios::failbit);
    return is;
}





} // namespace GEMSS

#endif // GEMSS_CONFIG_HPP
