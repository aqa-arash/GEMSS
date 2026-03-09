cmake_minimum_required(VERSION 3.15)
project(MultisphereLib VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Build Options ---
option(USE_OPENMP "Use OpenMP for parallelization" ON)
option(USE_MANIFOLD "Use Manifold for boolean operations" OFF)
option(USE_VTK "Enable VTK for visualization" OFF)

# Optimization settings
if(NOT CMAKE_BUILD_TYPE MATCHES Debug)
    if(MSVC)
        add_compile_options(/W4 /Ox /Oi /Ot)
    else()
        add_compile_options(-O3 -march=native)
    endif()
endif()

# --- Dependencies ---

# 1. Eigen (Required)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# 2. ZLIB (Required for cnpy)
find_package(ZLIB REQUIRED)

# 3. OpenMP (Optional)
if(USE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found and enabled.")
        add_definitions(-DHAVE_OPENMP)
    else()
        message(WARNING "OpenMP requested but not found. Compiling single-threaded.")
        set(USE_OPENMP OFF)
    endif()
endif()

# 4. VTK (Optional for Visualization)
if(USE_VTK)
    find_package(VTK COMPONENTS 
        CommonCore FiltersSources RenderingCore RenderingOpenGL2 InteractionStyle QUIET
    )
    if(VTK_FOUND)
        message(STATUS "VTK found: Visualization enabled.")
        add_definitions(-DHAVE_VTK)
    else()
        message(WARNING "VTK NOT found. Visualization disabled.")
        set(USE_VTK OFF)
    endif()
endif()

# 5. Manifold (Optional)
if(USE_MANIFOLD)
    find_package(manifold QUIET)
    if(manifold_FOUND)
        message(STATUS "Manifold found: Boolean operations enabled.")
        add_definitions(-DHAVE_MANIFOLD)
    else()
        message(WARNING "Manifold NOT found. Boolean logic disabled.")
    endif()
endif()

# --- Target Configuration ---

# A. Create a static library for cnpy (compiles cnpy.cpp)
add_library(cnpy_static STATIC include/cnpy.cpp)
target_include_directories(cnpy_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cnpy_static PRIVATE ZLIB::ZLIB)

# B. multisphere_lib (Header-only INTERFACE)
add_library(multisphere_lib INTERFACE)
target_include_directories(multisphere_lib INTERFACE 
                    ${CMAKE_CURRENT_SOURCE_DIR}
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/igl/
                    
)

# Link dependencies to the interface library
target_link_libraries(multisphere_lib INTERFACE 
    Eigen3::Eigen 
    cnpy_static
)

if(USE_OPENMP AND OpenMP_CXX_FOUND)
    target_link_libraries(multisphere_lib INTERFACE OpenMP::OpenMP_CXX)
endif()

if(VTK_FOUND)
    target_link_libraries(multisphere_lib INTERFACE ${VTK_LIBRARIES})
    vtk_module_autoinit(TARGETS multisphere_lib MODULES ${VTK_LIBRARIES})
endif()

if(manifold_FOUND AND USE_MANIFOLD)
    target_link_libraries(multisphere_lib INTERFACE manifold::manifold)
endif()

# --- Example App ---
add_executable(multisphere_app main.cpp)
add_executable(multisphere_mesh main_mesh.cpp)
add_executable(multisphere_collision main_collision.cpp)
add_executable(multisphere_boost main_boost.cpp)



target_link_libraries(multisphere_app PRIVATE multisphere_lib)
target_link_libraries(multisphere_mesh PRIVATE multisphere_lib)
target_link_libraries(multisphere_collision PRIVATE multisphere_lib)
target_link_libraries(multisphere_boost PRIVATE multisphere_lib)

