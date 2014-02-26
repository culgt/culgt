project(test_lattice)
cmake_minimum_required(VERSION 2.8)

include( ../commonsettings.cmake )

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



set(sources_cuda 
	test_LocalLink.cu
	test_GlobalLinkWithPatterns.cu
	../../lattice/cuda/CudaError.cu
	)

cuda_add_executable(test_lattice_cuda ../test_main.cpp ${sources_cuda} )
target_link_libraries(test_lattice_cuda pthread)
target_link_libraries(test_lattice_cuda gmock)
target_link_libraries(test_lattice_cuda gtest)

if( DEFINED allsources )
addWithDirectory( lattice "${sources_cuda}"  )
endif()



cuda_add_executable(test_lattice ../test_main.cpp ${sources} ${sources_cuda})
target_link_libraries(test_lattice pthread)
target_link_libraries(test_lattice gmock)
target_link_libraries(test_lattice gtest)
