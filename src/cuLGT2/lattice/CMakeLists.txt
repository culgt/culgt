project(test_lattice)
cmake_minimum_required(VERSION 2.8)

include( ../addWithDirectory.cmake )

include_directories($ENV{GMOCK_HOME}/include $ENV{GMOCK_HOME}/gtest/include)
link_directories($ENV{GMOCK_HOME}/mybuild $ENV{GMOCK_HOME}/gtest/mybuild)
add_definitions()
set(CMAKE_CXX_FLAGS "${CMAXE_CXX_FLAGS} -Wall")

set(sources 
	test_LocalLink.cc
	su3/ParameterizationMediatorSU3_Real12_Real18_test.cc 
	su3/SU3.cc
	configuration_patterns/test_StandardPattern.cc
	configuration_patterns/test_GPUPattern.cc
	configuration_patterns/test_GPUPatternParityPriority.cc
	test_GlobalLink.cc 
	test_GlobalLinkWithPatterns.cc
	test_LocalGlobalLinkCompatibility.cc
	)
add_executable(test_lattice_host ../test_main.cpp ${sources})
target_link_libraries(test_lattice_host pthread)
target_link_libraries(test_lattice_host gmock)
target_link_libraries(test_lattice_host gtest)


if( DEFINED allsources )
addWithDirectory( lattice "${sources}" )
endif()

find_package(CUDA)
include(FindCUDA)

set(sources_cuda 
	test_LocalLink.cu
	test_GlobalLinkWithPatterns.cu
	../../lattice/cuda/CudaError.cu
	)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-arch=sm_20;)

CUDA_ADD_EXECUTABLE(test_lattice_cuda ../test_main.cpp ${sources_cuda} )
target_link_libraries(test_lattice_cuda pthread)
target_link_libraries(test_lattice_cuda gmock)
target_link_libraries(test_lattice_cuda gtest)

if( DEFINED allsources )
addWithDirectory( lattice "${sources_cuda}"  )
endif()

CUDA_ADD_EXECUTABLE(test_lattice ../test_main.cpp ${sources} ${sources_cuda})
target_link_libraries(test_lattice pthread)
target_link_libraries(test_lattice gmock)
target_link_libraries(test_lattice gtest)
