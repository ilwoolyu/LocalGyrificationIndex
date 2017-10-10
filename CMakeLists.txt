cmake_minimum_required(VERSION 2.8)

project(Gyrification)

# CMAKE
include_directories(${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})

# GenerateCLP
find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})
find_package(GenerateCLP REQUIRED)
include(${GenerateCLP_USE_FILE})
generateCLP(WRAPPER_XML wrapper/wrapper.xml)
set(WRAPPER_MAIN wrapper/main.cpp)

# MESH library
find_package(MeshLib REQUIRED)
include_directories(${Mesh_INCLUDE_DIRS} ${Mesh_INCLUDE_DIRS}/Util)
link_directories(${Mesh_LIBRARY_DIRS})

# sources
include_directories(src)
set(SRC_FILES src/Gyrification.cpp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ~/bin)

add_executable(Gyrification ${SRC_FILES} ${WRAPPER_XML} ${WRAPPER_MAIN})
target_link_libraries(Gyrification Mesh)
