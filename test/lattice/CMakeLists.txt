project(test_lattice)

cmake_minimum_required(VERSION 2.8)

include( ../commonsettings.cmake )



set(sources 
	test_LocalLink.cc
	test_GlobalLink.cc 
	test_GlobalLinkWithPatterns.cc
	test_LocalGlobalLinkCompatibility.cc
	test_LinkToGaugefieldConverter.cc
	parameterization_types/test_ParameterizationMediatorSU3_Real12_Real18.cc 
	parameterization_types/test_ParameterizationMediatorSU2_Real4_Real8.cc 
	parameterization_types/test_ParameterizationMediatorSU3_Vector4_Complex9.cc 
	parameterization_types/test_SUNRealFull.cc
	parameterization_types/test_SU3Real12.cc
	configuration_patterns/test_StandardPattern.cc
	configuration_patterns/test_GPUPattern.cc
	configuration_patterns/test_GPUPatternTimeslice.cc
	configuration_patterns/test_GPUPatternTimesliceParityPriority.cc
	configuration_patterns/test_GPUPatternParityPriority.cc
	filetypes/test_limereader.cc
	filetypes/test_LinkFileVogt.cc
	filetypes/test_LinkFileVogtWithPattern.cc
	filetypes/test_LinkFile.cc
	filetypes/test_LinkFileManager.cc
	site_indexing/test_SiteIndex.cc
	site_indexing/test_SiteCoord.cc
	site_indexing/test_SiteCompatibilityWithConfigurationPatterns.cc
	site_indexing/test_SiteIndexSiteCoordCompatibility.cc
	)
	
add_library( lib_lattice_host ${sources} )
target_link_libraries( lib_lattice_host pthread )
target_link_libraries( lib_lattice_host gmock )
target_link_libraries( lib_lattice_host gtest )
target_link_libraries( lib_lattice_host liblime.a )
target_link_libraries( lib_lattice_host tinyxml )

add_executable( test_lattice_host ../test_main.cpp )
target_link_libraries_whole_archive( test_lattice_host lib_lattice_host )



set(sources_cuda 
	test_LocalLink.cu
	site_indexing/test_SiteIndex.cu
	test_LocalLinkWithSU2Vector4.cu
	test_GlobalLinkWithPatterns.cu
	test_GaugeConfiguration.cu
	test_GaugeConfigurationWithPattern.cu
	test_LatticeDimension.cu
	filetypes/test_LinkFileCompatibility.cu
	test_KernelSetup.cu
	parameterization_types/test_SU2Vector4.cu
	site_indexing/test_SiteNeighbourTableManager.cu
	)

cuda_add_library( lib_lattice_cuda ${sources_cuda} )
target_link_libraries( lib_lattice_cuda culgt_cudacommon )
target_link_libraries( lib_lattice_cuda pthread )
target_link_libraries( lib_lattice_cuda gmock )
target_link_libraries( lib_lattice_cuda gtest )
target_link_libraries( lib_lattice_cuda liblime.a )
target_link_libraries( lib_lattice_cuda tinyxml )

cuda_add_executable( test_lattice_cuda ../test_main.cpp )
target_link_libraries_whole_archive( test_lattice_cuda lib_lattice_cuda )



cuda_add_executable(test_lattice ../test_main.cpp )
target_link_libraries_whole_archive( test_lattice lib_lattice_cuda )
target_link_libraries_whole_archive( test_lattice lib_lattice_host )

