project(test_lattice)

cmake_minimum_required(VERSION 2.8)

include( ../commonsettings.cmake )

#set( CMAKE_VERBOSE_MAKEFILE on )

set(sources 
	test_LocalLink.cc
	parameterization_types/ParameterizationMediatorSU3_Real12_Real18_test.cc 
	parameterization_types/test_ParameterizationMediatorSU2_Real4_Real8.cc 
	parameterization_types/SU3.cc
	parameterization_types/test_SUNRealFull.cc
	parameterization_types/test_SU3Real12.cc
	configuration_patterns/test_StandardPattern.cc
	configuration_patterns/test_GPUPattern.cc
	configuration_patterns/test_GPUPatternTimeslice.cc
	configuration_patterns/test_GPUPatternTimesliceParityPriority.cc
	configuration_patterns/test_GPUPatternParityPriority.cc
	test_GlobalLink.cc 
	test_GlobalLinkWithPatterns.cc
	test_LocalGlobalLinkCompatibility.cc
	test_LinkFileVogt.cc
	test_LinkFileVogtWithPattern.cc
	test_LinkFile.cc
	test_KernelSetup.cc
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
	test_LocalLinkWithSU2Vector4.cu
	test_GlobalLinkWithPatterns.cu
	../../lattice/cuda/CudaError.cu
	test_GaugeConfiguration.cc
	test_GaugeConfigurationWithPattern.cu
	test_LatticeDimension.cu
	../cudacommon/test_DeviceCommunicator.cu
	parameterization_types/test_SU2Vector4.cu
	)

cuda_add_executable(test_lattice_cuda ../test_main.cpp ${sources_cuda} )
target_link_libraries(test_lattice_cuda pthread)
target_link_libraries(test_lattice_cuda gmock)
target_link_libraries(test_lattice_cuda gtest)
target_link_libraries(test_lattice_cuda cudart)

if( DEFINED allsources )
addWithDirectory( lattice "${sources_cuda}"  )
endif()



cuda_add_executable(test_lattice ../test_main.cpp ${sources} ${sources_cuda})
target_link_libraries(test_lattice pthread)
target_link_libraries(test_lattice gmock)
target_link_libraries(test_lattice gtest)





set(sourcesSpecial 
	test_GlobalLinkWithPatterns.cc
	)
add_executable(test_lattice_special ../test_main.cpp ${sourcesSpecial})
target_link_libraries(test_lattice_special pthread)
target_link_libraries(test_lattice_special gmock)
target_link_libraries(test_lattice_special gtest)