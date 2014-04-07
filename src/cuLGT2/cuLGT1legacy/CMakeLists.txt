project(test_cuLGT1legacy)
cmake_minimum_required(VERSION 2.8)

include( ../commonsettings.cmake )

set(sources 
	test_SiteIndex.cc
	test_SiteCoord.cc
	test_SiteCompatibilityWithConfigurationPatterns.cc
	test_SiteIndexSiteCoordCompatibility.cc
	)
add_executable(test_cuLGT1legacy_host ../test_main.cpp  ${sources})
target_link_libraries(test_cuLGT1legacy_host pthread)
target_link_libraries(test_cuLGT1legacy_host gmock)
target_link_libraries(test_cuLGT1legacy_host gtest)

if( DEFINED allsources )
addWithDirectory( cuLGT1legacy "${sources}" )
endif()

set(sources_cuda 
	test_SiteNeighbourTableManager.cu
	)

cuda_add_executable(test_cuLGT1legacy_cuda ../test_main.cpp ${sources_cuda} )
target_link_libraries(test_cuLGT1legacy_cuda pthread)
target_link_libraries(test_cuLGT1legacy_cuda gmock)
target_link_libraries(test_cuLGT1legacy_cuda gtest)
target_link_libraries(test_cuLGT1legacy_cuda cudart)

if( DEFINED allsources )
addWithDirectory( cuLGT1legacy "${sources_cuda}"  )
endif()




cuda_add_executable(test_cuLGT1legacy ../test_main.cpp ${sources} ${sources_cuda})
target_link_libraries(test_cuLGT1legacy pthread)
target_link_libraries(test_cuLGT1legacy gmock)
target_link_libraries(test_cuLGT1legacy gtest)