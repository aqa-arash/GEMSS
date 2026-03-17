cmake_minimum_required(VERSION 3.15)
project(MultisphereLib VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# --- Build Options ---
option(USE_OPENMP "Use OpenMP for parallelization" ON)
option(MULTISPHERE_DEBUG "Enable console debug prints (MULTISPHERE_DEBUG macro)" OFF)
option(MULTISPHERE_BUILD_EXAMPLES "Build example executables" ON)

# --- Dependencies ---
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

if(USE_OPENMP)
    find_package(OpenMP)
    if(NOT OpenMP_CXX_FOUND)
        message(WARNING "OpenMP requested but not found. Compiling single-threaded.")
        set(USE_OPENMP OFF)
    endif()
endif()

# --- Library Target ---
# Multisphere Core (Pure INTERFACE / Header-Only)
add_library(multisphere_lib INTERFACE)

target_include_directories(multisphere_lib INTERFACE 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/thirdparty>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/thirdparty/igl>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(multisphere_lib INTERFACE Eigen3::Eigen)

if(USE_OPENMP AND OpenMP_CXX_FOUND)
    target_link_libraries(multisphere_lib INTERFACE OpenMP::OpenMP_CXX)
    target_compile_definitions(multisphere_lib INTERFACE HAVE_OPENMP)
endif()

# Toggle console prints based on the specific option
if(MULTISPHERE_DEBUG)
    target_compile_definitions(multisphere_lib INTERFACE MULTISPHERE_DEBUG)
endif()

# --- Example Applications ---
if(MULTISPHERE_BUILD_EXAMPLES)
    # UNCONDITIONAL optimization flags for example executables
    if(MSVC)
        set(OPT_FLAGS /Ox /Oi /Ot)
    else()
        set(OPT_FLAGS -O3 -march=native)
    endif()

    macro(add_multisphere_example target_name source_file)
        if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${source_file}")
            add_executable(${target_name} ${source_file})
            target_link_libraries(${target_name} PRIVATE multisphere_lib)
            # Force optimizations regardless of CMake build type
            target_compile_options(${target_name} PRIVATE ${OPT_FLAGS})
        endif()
    endmacro()

    add_multisphere_example(multisphere_app main.cpp)
    add_multisphere_example(multisphere_mesh main_mesh.cpp)
    add_multisphere_example(multisphere_collision main_collision.cpp)
    add_multisphere_example(multisphere_boost main_boost.cpp)
endif()

# --- Installation Rules ---
include(GNUInstallDirs)
install(TARGETS multisphere_lib
    EXPORT MultisphereTargets
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})