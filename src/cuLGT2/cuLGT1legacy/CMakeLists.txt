project(test_cuLGT1legacy)
cmake_minimum_required(VERSION 2.8)

include( ../commonsettings.cmake )

set(sources 
	test_SiteIndex.cc
	test_SiteCoord.cc
	test_SiteCompatibilityWithConfigurationPatterns.cc
	)
add_executable(test_cuLGT1legacy ../test_main.cpp  ${sources})
target_link_libraries(test_cuLGT1legacy pthread)
target_link_libraries(test_cuLGT1legacy gmock)
target_link_libraries(test_cuLGT1legacy gtest)

if( DEFINED allsources )
addWithDirectory( cuLGT1legacy "${sources}" )
endif()

