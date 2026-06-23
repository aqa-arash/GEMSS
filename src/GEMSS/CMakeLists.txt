cmake_minimum_required(VERSION 3.15)
project(GEMSS LANGUAGES CXX)

# Needed to populate CMAKE_INSTALL_INCLUDEDIR
include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# --- Build Options ---
option(GEMSS_USE_OPENMP "Use OpenMP for parallelization" ON)
option(GEMSS_DEBUG "Enable console debug prints (GEMSS_DEBUG macro)" OFF)

# --- Default Build Type to Release ---
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "No build type selected, defaulting to Release for maximum performance.")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()

# --- Dependencies ---
find_package(Eigen3 REQUIRED NO_MODULE)

if(GEMSS_USE_OPENMP)
    find_package(OpenMP)
    if(NOT OpenMP_CXX_FOUND)
        message(WARNING "OpenMP requested but not found. Compiling single-threaded.")
        set(GEMSS_USE_OPENMP OFF)
    endif()
endif()

# --- Library Target ---

add_library(GEMSS INTERFACE)
# Create alias
add_library(GEMSS::GEMSS ALIAS GEMSS)


target_include_directories(GEMSS INTERFACE 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/thirdparty>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/igl>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/GEMSS>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/GEMSS/thirdparty>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/GEMSS/thirdparty/igl>
)
target_link_libraries(GEMSS INTERFACE Eigen3::Eigen)

if(GEMSS_USE_OPENMP AND OpenMP_CXX_FOUND)
    target_link_libraries(GEMSS INTERFACE OpenMP::OpenMP_CXX)
    target_compile_definitions(GEMSS INTERFACE HAVE_OPENMP)
endif()

if(GEMSS_DEBUG)
    target_compile_definitions(GEMSS INTERFACE GEMSS_DEBUG)
endif()

# --- Installation Rules if GEMSS is the main project being built
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CMakePackageConfigHelpers)

    # Install the library target and associate it with the GEMSSTargets export set
    install(TARGETS GEMSS
        EXPORT GEMSSTargets
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/GEMSS
    )


    # Install header files (Avoids installing CMakeLists.txt or binaries into the include folder)
    install(DIRECTORY .
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/GEMSS
        FILES_MATCHING 
        PATTERN "*.h" 
        PATTERN "*.hpp" 
        PATTERN "*.cpp"
        PATTERN "*.c"
    )

    # 1. Export the targets to create GEMSSTargets.cmake
    install(EXPORT GEMSSTargets
        FILE GEMSSTargets.cmake
        NAMESPACE GEMSS::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GEMSS
    )

    # 2. Configure the GEMSSConfig.cmake file from your template
    configure_package_config_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/GEMSSConfig.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/GEMSSConfig.cmake"
        INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/GEMSS"
    )

    # 3. Create a version file (Good practice for package management)
    set(GEMSS_VERSION "1.0.0")
    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/GEMSSConfigVersion.cmake"
        VERSION ${GEMSS_VERSION}
        COMPATIBILITY AnyNewerVersion
    )

    # 4. Install the generated configuration files
    install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/GEMSSConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/GEMSSConfigVersion.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/GEMSS"
    )
endif()