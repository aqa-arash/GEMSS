cmake_minimum_required(VERSION 3.15)

project(GEMSSExamples LANGUAGES CXX)

# Option to build GUI
option(BUILD_GUI "Build GEMSS Qt/VTK GUI" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Add GEMSS header-only library
add_subdirectory(GEMSS)


# Gather all .cpp files except gemss_gui.cpp if BUILD_GUI is OFF
file(GLOB cpp_files "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
set(GUI_SRC "${CMAKE_CURRENT_SOURCE_DIR}/gemss_gui.cpp")
set(FILTERED_CPP_FILES "")
foreach(src_file ${cpp_files})
    if(NOT src_file STREQUAL "${GUI_SRC}")
        list(APPEND FILTERED_CPP_FILES ${src_file})
    endif()
endforeach()

# Build regular executables
foreach(src_file ${FILTERED_CPP_FILES})
    get_filename_component(target_name ${src_file} NAME_WE)
    add_executable(${target_name} ${src_file})
    target_link_libraries(${target_name} PRIVATE GEMSS)
    target_compile_options(${target_name} PRIVATE -O3 -march=native)
endforeach()

# Build GUI if requested
if(BUILD_GUI)
    # Find Qt and VTK (VTK 9+ modules)
    find_package(Qt5 COMPONENTS Widgets REQUIRED)
    find_package(VTK COMPONENTS 
        CommonCore CommonDataModel FiltersSources FiltersCore FiltersGeneral 
        InteractionStyle RenderingCore RenderingOpenGL2 RenderingContextOpenGL2 
        RenderingFreeType RenderingGL2PSOpenGL2 GUISupportQt RenderingQt ViewsQt 
        IOGeometry IOImport IOPLY RenderingAnnotation
    )
    if (NOT VTK_FOUND)
        message(FATAL_ERROR "VTK not found!")
    endif()
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTOUIC ON)
    set(CMAKE_AUTORCC ON)
    add_executable(gemss_gui gemss_gui.cpp)
    target_link_libraries(gemss_gui PRIVATE GEMSS Qt5::Widgets ${VTK_LIBRARIES})
    target_compile_options(gemss_gui PRIVATE -O3 -march=native)
    if (VTK_VERSION VERSION_GREATER_EQUAL "9.0")
        include(${VTK_USE_FILE})
        vtk_module_autoinit(
            TARGETS gemss_gui
            MODULES ${VTK_LIBRARIES}
        )
    endif()
endif()




