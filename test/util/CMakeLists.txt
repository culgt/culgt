project(test_util)

cmake_minimum_required(VERSION 2.8)

include_directories(../../include)
include( ../commonsettings.cmake )

#set( CMAKE_VERBOSE_MAKEFILE on )

set(sources 
	performance/test_AutotuneManager.cc
	)
add_executable(test_util_host ../test_main.cpp ${sources})
target_link_libraries(test_util_host pthread)
target_link_libraries(test_util_host gmock)
target_link_libraries(test_util_host gtest)

if( DEFINED allsources )
addWithDirectory( util "${sources}" )
endif()



set(sources_cuda
    reduction/test_Reduction.cu
	)

cuda_add_executable(test_util_cuda ../test_main.cpp ${sources_cuda} )
target_link_libraries(test_util_cuda pthread)
target_link_libraries(test_util_cuda gmock)
target_link_libraries(test_util_cuda gtest)
target_link_libraries(test_util_cuda cudart)

if( DEFINED allsources )
addWithDirectory( util "${sources_cuda}"  )
endif()



cuda_add_executable(test_util ../test_main.cpp ${sources} ${sources_cuda})
target_link_libraries(test_util pthread)
target_link_libraries(test_util gmock)
target_link_libraries(test_util gtest)



