cmake_minimum_required(VERSION 3.20)
project(ToyInterpolation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# Find packages
find_package(ROOT REQUIRED)
find_package(Eigen3 REQUIRED)
# find_package(Boost REQUIRED COMPONENTS program_options)
find_package(Boost REQUIRED COMPONENTS program_options python)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Compile every .cxx file in src/ into object files
file(GLOB SRC_FILES ${CMAKE_SOURCE_DIR}/src/*.cxx)
add_library(toy_interpolation_objects OBJECT ${SRC_FILES})
target_link_libraries(toy_interpolation_objects PUBLIC Eigen3::Eigen) 
target_include_directories(toy_interpolation_objects PUBLIC ${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS})

# Add executables
add_executable(interpolation interpolation.cxx)
add_executable(honda_flux_handling honda_flux_handling.cxx)
add_executable(honda_flux_handling3D honda_flux_handling3D.cxx)
add_executable(honda_flux_handling3D_spline honda_flux_handling3D_spline.cxx)

# Link libraries
target_link_libraries(interpolation PRIVATE ROOT::Core ROOT::Hist ROOT::Gpad Eigen3::Eigen)
target_link_libraries(honda_flux_handling PRIVATE Boost::program_options toy_interpolation_objects ROOT::Core ROOT::Hist ROOT::Gpad Eigen3::Eigen)
target_link_libraries(honda_flux_handling3D PRIVATE Boost::program_options toy_interpolation_objects ROOT::Core ROOT::Hist ROOT::Gpad Eigen3::Eigen)
target_link_libraries(honda_flux_handling3D_spline PRIVATE Boost::program_options toy_interpolation_objects ROOT::Core ROOT::Hist ROOT::Gpad Eigen3::Eigen)

# Add shared library for Python interface
add_library(hkkm_interpolation MODULE py_interface.cxx)

# Link Boost.Python to the shared library
target_link_libraries(hkkm_interpolation PRIVATE Boost::python Python3::Python toy_interpolation_objects Eigen3::Eigen ROOT::Core ROOT::Hist)

# Set the output properties for the shared library
set_target_properties(hkkm_interpolation PROPERTIES PREFIX "" SUFFIX ".so")
