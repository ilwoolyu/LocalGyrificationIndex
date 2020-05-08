set(CLI11_VER "1.9.0" CACHE STRING "CLI11 Version: see https://github.com/CLIUtils/CLI11/releases")
mark_as_advanced(CLI11_VER)
set(CLI11_DIRS ${CMAKE_BINARY_DIR}/CMakeFiles/CLI11)
file(DOWNLOAD https://github.com/CLIUtils/CLI11/releases/download/v${CLI11_VER}/CLI11.hpp ${CLI11_DIRS}/CLI11.hpp)
include_directories(${CLI11_DIRS})

include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(
	Gyrification
	${CMAKE_SOURCE_DIR}/src/Gyrification.cpp
	${CMAKE_SOURCE_DIR}/wrapper/main.cpp
)

if (OPENMP_FOUND)
	target_compile_definitions(Gyrification PRIVATE _USE_OPENMP)
endif()
target_compile_definitions(Gyrification PRIVATE LGI_VERSION=\"v${Gyrification_VERSION}\")

target_link_libraries(Gyrification Mesh)
