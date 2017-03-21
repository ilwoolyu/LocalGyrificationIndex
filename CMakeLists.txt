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
set(MESH_LIB_DIR /NIRAL/work/ilwoolyu/MeshLib)
include_directories(${MESH_LIB_DIR}/include)
link_directories(${MESH_LIB_DIR}/lib)

# sources
include_directories(src)
set(SRC_FILES src/Gyrification.cpp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

add_executable(Gyrification ${SRC_FILES} ${WRAPPER_XML} ${WRAPPER_MAIN})
#add_executable(Gyrification ${SRC_FILES} ${WRAPPER_MAIN})
target_link_libraries(Gyrification Mesh)
