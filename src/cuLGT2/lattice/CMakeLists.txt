project(GaugeConfiguration)
cmake_minimum_required(VERSION 2.8)

include_directories($ENV{GMOCK_HOME}/include $ENV{GMOCK_HOME}/gtest/include)
link_directories($ENV{GMOCK_HOME}/mybuild $ENV{GMOCK_HOME}/gtest/mybuild)
add_definitions()
set(CMAKE_CXX_FLAGS "${CMAXE_CXX_FLAGS} -Wall")

set(sources 
	../test_main.cpp 
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
add_executable(test_host ${sources})
target_link_libraries(test_host pthread)
target_link_libraries(test_host gmock)
target_link_libraries(test_host gtest)


find_package(CUDA)
include(FindCUDA)

set(sources_cuda 
	../test_main.cpp 
	test_LocalLink.cu
	test_GlobalLinkWithPatterns.cu
	../../lattice/cuda/CudaError.cu
	)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-arch=sm_20;)

CUDA_ADD_EXECUTABLE(test_cuda ${sources_cuda} )
target_link_libraries(test_cuda pthread)
target_link_libraries(test_cuda gmock)
target_link_libraries(test_cuda gtest)


CUDA_ADD_EXECUTABLE(test ${sources} ${sources_cuda})
target_link_libraries(test pthread)
target_link_libraries(test gmock)
target_link_libraries(test gtest)
