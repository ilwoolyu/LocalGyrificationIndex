cmake_minimum_required(VERSION 3.5)

project(Gyrification)

# Set a default build type if none was specified
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        message(STATUS "Setting build type to 'Release' as none was specified.")
        set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
        mark_as_advanced(CMAKE_BUILD_TYPE)
        # Set the possible values of build type for cmake-gui
        set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# GenerateCLP
find_package(SlicerExecutionModel REQUIRED)
include(${GenerateCLP_USE_FILE})
set(WRAPPER_MAIN wrapper/main.cpp)
generateCLP(WRAPPER_MAIN ${CMAKE_SOURCE_DIR}/wrapper/wrapper.xml)

# MESH library
find_package(MeshLib REQUIRED)

# sources
include_directories(src)
set(SRC_FILES ${CMAKE_SOURCE_DIR}/src/Gyrification.cpp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ~/bin)

add_executable(Gyrification ${SRC_FILES} ${WRAPPER_MAIN})
target_link_libraries(Gyrification Mesh)
