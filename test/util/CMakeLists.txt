project(test_util)
cmake_minimum_required(VERSION 2.8)

include( ../commonsettings.cmake )

set(sources 
	performance/test_AutotuneManager.cc
	)
add_library(lib_util_host ${sources})
target_link_libraries(lib_util_host pthread)
target_link_libraries(lib_util_host gmock)
target_link_libraries(lib_util_host gtest)

add_executable( test_util_host ../test_main.cpp )
target_link_libraries_whole_archive( test_util_host lib_util_host )



set(sources_cuda
    reduction/test_Reduction.cu
	)

cuda_add_library(lib_util_cuda  ${sources_cuda} )
target_link_libraries(lib_util_cuda pthread)
target_link_libraries(lib_util_cuda gmock)
target_link_libraries(lib_util_cuda gtest)

cuda_add_executable( test_util_cuda ../test_main.cpp )
target_link_libraries_whole_archive( test_util_cuda lib_util_cuda )



cuda_add_executable(test_util ../test_main.cpp )
target_link_libraries_whole_archive( test_util lib_util_host )
target_link_libraries_whole_archive( test_util lib_util_cuda )
