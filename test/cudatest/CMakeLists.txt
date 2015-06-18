# IMPORTANT: This demonstrates the use of the cuda_test code. For demonstration some tests fail. Therefore
# we do not include this subdirectory in the full test.

project(cuda_test_helper_test)
cmake_minimum_required(VERSION 2.8)



include( ../commonsettings.cmake )
 
set(sources 
	test_example_gtest_plugin_host.cc
	)
add_library(lib_cuda_test_helper_host ${sources} )
target_link_libraries(lib_cuda_test_helper_host pthread)
target_link_libraries(lib_cuda_test_helper_host gmock)
target_link_libraries(lib_cuda_test_helper_host gtest)
 
add_executable( test_cuda_test_helper_host ../test_main.cpp )
target_link_libraries_whole_archive( test_cuda_test_helper_host lib_cuda_test_helper_host )



set(sources_cuda 
	test_example_gtest_plugin.cu
	)

cuda_add_library(lib_cuda_test_helper_cuda ${sources_cuda} )
target_link_libraries(lib_cuda_test_helper_cuda pthread)
target_link_libraries(lib_cuda_test_helper_cuda gmock)
target_link_libraries(lib_cuda_test_helper_cuda gtest)

add_executable( test_cuda_test_helper_cuda ../test_main.cpp )
target_link_libraries_whole_archive( test_cuda_test_helper_cuda lib_cuda_test_helper_cuda )



add_executable( test_cuda_test_helper ../test_main.cpp )
target_link_libraries_whole_archive( test_cuda_test_helper lib_cuda_test_helper_host )
target_link_libraries_whole_archive( test_cuda_test_helper lib_cuda_test_helper_cuda )