cmake_minimum_required(VERSION 3.20)
project(ToyInterpolation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find packages
find_package(ROOT REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Compile every .cxx file in src/ into object files
file(GLOB SRC_FILES ${CMAKE_SOURCE_DIR}/src/*.cxx)
add_library(toy_interpolation_objects OBJECT ${SRC_FILES})
# link_libraries(toy_interpolation_objects PUBLIC ROOT::Core ROOT::Hist ROOT::Gpad) 
target_include_directories(toy_interpolation_objects PUBLIC ${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS})

# Add executables
add_executable(interpolation interpolation.cxx)
add_executable(honda_flux_handling honda_flux_handling.cxx)
add_executable(honda_flux_handling3D honda_flux_handling3D.cxx)

# Link libraries
target_link_libraries(interpolation PRIVATE ROOT::Core ROOT::Hist ROOT::Gpad Eigen3::Eigen)
target_link_libraries(honda_flux_handling PRIVATE Boost::program_options toy_interpolation_objects ROOT::Core ROOT::Hist ROOT::Gpad Eigen3::Eigen)
target_link_libraries(honda_flux_handling3D PRIVATE Boost::program_options toy_interpolation_objects ROOT::Core ROOT::Hist ROOT::Gpad Eigen3::Eigen)