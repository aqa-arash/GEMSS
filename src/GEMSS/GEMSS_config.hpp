#ifndef GEMSS_CONFIG_HPP
#define GEMSS_CONFIG_HPP


#include <Eigen/Core>
#include <istream>
#include <sstream>
#include <limits>
#include <ostream>
#include <string>

namespace GEMSS {


struct MultisphereConfig {
    // --- Required for mesh reconstruction ---
    int div = 100;               // grid resolution
    bool volumetric_div = false; // voxel size: cbrt(volume) / div, at least 2 voxels across the smallest extent (true); smallest extent / div (false)
    int padding = 2;             // grid padding
    float minimum_radius_real = 0.0f; // min sphere radius (real units)
    bool confine_mesh = false;   // confine to mesh boundary

    // --- Peak detection and filtering ---
    float min_center_distance_rel = 1.4f; // min center distance (relative) computed as k <= sqrt(min_radius)
    int search_window = 2;       // neighbor distance for peaks
    float radius_offset_vox = 0.5f; // radius correction
    int min_radius_vox = 2;      // min sphere radius (voxels)
    float medial_flux_threshold = 0.05f; // peaks are searched where the average outward flux of grad(distance) < -threshold (medial axis); < 0 searches every voxel


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
       //<< " conserve_mass " << c.conserve_mass
       //<< " parent_representation_offset " << c.parent_representation_offset
       //<< " fracture_plane_offset " << c.fracture_plane_offset
       << " end";
    os.precision(precision);
    return os;
}

// Reads "name value" pairs until "end". Each value is matched to its setting by name, so the order does not
// matter, a setting missing from the input keeps its default, and an unknown setting (e.g. from another GEMSS
// version) is skipped with its value.
inline std::istream& operator>>(std::istream& is, MultisphereConfig& c)
{
    std::string word;
    if (!(is >> word) || word != "MultisphereConfig") { is.setstate(std::ios::failbit); return is; }
    std::string key, token;
    while (is >> key && key != "end")
    {
        if (!(is >> token)) break;
        std::istringstream value(token);
        if (key == "div") value >> c.div;
        else if (key == "volumetric_div") value >> c.volumetric_div;
        else if (key == "padding") value >> c.padding;
        else if (key == "minimum_radius_real") value >> c.minimum_radius_real;
        else if (key == "confine_mesh") value >> c.confine_mesh;
        else if (key == "min_center_distance_rel") value >> c.min_center_distance_rel;
        else if (key == "search_window") value >> c.search_window;
        else if (key == "radius_offset_vox") value >> c.radius_offset_vox;
        else if (key == "min_radius_vox") value >> c.min_radius_vox;
        else if (key == "medial_flux_threshold") value >> c.medial_flux_threshold;
        else if (key == "precision_target") value >> c.precision_target;
        else if (key == "max_spheres") value >> c.max_spheres;
        else if (key == "persistence") value >> c.persistence;
        else if (key == "compute_physics") value >> c.compute_physics;
        else if (key == "is_2d") value >> c.is_2d;
        else if (key == "density") value >> c.density;
        else if (key == "prune_isolated_spheres") value >> c.prune_isolated_spheres;
        else if (key == "show_progress") value >> c.show_progress;
    }
    if (key != "end") is.setstate(std::ios::failbit);
    return is;
}





} // namespace GEMSS

#endif // GEMSS_CONFIG_HPP
