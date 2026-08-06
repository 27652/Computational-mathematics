cmake_minimum_required(VERSION 3.18)

project(ROPTLIB_TestRPNCG LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_library(roptlib_test_rpncg_obj OBJECT)

target_sources(roptlib_test_rpncg_obj PRIVATE
	Manifolds/Element.cpp
	Manifolds/Manifold.cpp
	Manifolds/Stiefel.cpp
	Others/BlasLapackCppWrapper.cpp
	Others/ForDebug.cpp
	Others/randgen.cpp
	Others/Timer.cpp
	Problems/StieSPCA.cpp
	Solvers/Solvers.cpp
	Solvers/SolversNSM.cpp
	Solvers/RPNCG.cpp
	test/mainRPNCG.cpp
	test/testRPNCG.cpp

Manifolds/SmartSpace.cpp
Others/SparseMatrix.cpp
Others/SparseBLAS/nist_spblas.cpp
Problems/Problem.cpp

)

target_include_directories(roptlib_test_rpncg_obj PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}/Manifolds
	${CMAKE_CURRENT_SOURCE_DIR}/Others
	${CMAKE_CURRENT_SOURCE_DIR}/Problems
	${CMAKE_CURRENT_SOURCE_DIR}/Solvers
	${CMAKE_CURRENT_SOURCE_DIR}/test
	${CMAKE_CURRENT_SOURCE_DIR}/cwrapper
	${CMAKE_CURRENT_SOURCE_DIR}/cwrapper/blas
	${CMAKE_CURRENT_SOURCE_DIR}/cwrapper/lapack
)

if(MSVC)
	target_compile_options(roptlib_test_rpncg_obj PRIVATE /permissive-)
endif()


add_executable(testRPNCG
    $<TARGET_OBJECTS:roptlib_test_rpncg_obj>
)
target_compile_options(roptlib_test_rpncg_obj PRIVATE
    -include cstring
)

set_target_properties(testRPNCG PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)


target_link_libraries(testRPNCG PRIVATE
    LAPACK::LAPACK
    BLAS::BLAS
)

