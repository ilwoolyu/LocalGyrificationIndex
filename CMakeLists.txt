cmake_minimum_required(VERSION 3.5)

project(Gyrification)

set (CMAKE_CXX_STANDARD 11)

set(Gyrification_VERSION_MAJOR 1)
set(Gyrification_VERSION_MINOR 2)
set(Gyrification_VERSION_PATCH 2)
set(Gyrification_VERSION
${Gyrification_VERSION_MAJOR}.${Gyrification_VERSION_MINOR}.${Gyrification_VERSION_PATCH})

# Set a default build type if none was specified
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        message(STATUS "Setting build type to 'Release' as none was specified.")
        set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
        mark_as_advanced(CMAKE_BUILD_TYPE)
        # Set the possible values of build type for cmake-gui
        set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

find_package(MeshLib REQUIRED)
find_package(OpenMP)

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# add subdiretories
add_subdirectory(wrapper)
