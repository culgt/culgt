project(test_lattice)

cmake_minimum_required(VERSION 2.8)

include_directories(../../include)
include( ../commonsettings.cmake )

#set( CMAKE_VERBOSE_MAKEFILE on )

set(sources 
	)
add_executable(test_gaugefixing_host ../test_main.cpp ${sources})
target_link_libraries(test_gaugefixing_host pthread)
target_link_libraries(test_gaugefixing_host gmock)
target_link_libraries(test_gaugefixing_host gtest)

if( DEFINED allsources )
addWithDirectory( gaugefixing "${sources}" )
endif()



set(sources_cuda 
	)

cuda_add_executable(test_gaugefixing_cuda ../test_main.cpp ${sources_cuda} )
target_link_libraries(test_gaugefixing_cuda pthread)
target_link_libraries(test_gaugefixing_cuda gmock)
target_link_libraries(test_gaugefixing_cuda gtest)
target_link_libraries(test_gaugefixing_cuda cudart)

if( DEFINED allsources )
addWithDirectory( gaugefixing "${sources_cuda}"  )
endif()



cuda_add_executable(test_gaugefixing ../test_main.cpp ${sources} ${sources_cuda})
target_link_libraries(test_gaugefixing pthread)
target_link_libraries(test_gaugefixing gmock)
target_link_libraries(test_gaugefixing gtest)



