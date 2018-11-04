include(${GenerateCLP_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR}/src)

generateCLP(Gyrification_WRAPPER ${CMAKE_SOURCE_DIR}/wrapper/wrapper.xml)

add_executable(
	Gyrification
	${CMAKE_SOURCE_DIR}/src/Gyrification.cpp
	${CMAKE_SOURCE_DIR}/wrapper/main.cpp ${Gyrification_WRAPPER}
)

if (OPENMP_FOUND)
	target_compile_definitions(Gyrification PRIVATE _USE_OPENMP)
endif()

target_link_libraries(Gyrification Mesh)
