project(test_eigensolver)
cmake_minimum_required(VERSION 2.8)

include( ../commonsettings.cmake )



set(sources 
	test_PowerIterator.cc
	test_Cardano3x3.cc
	test_Matrix.cc
	)
add_library(lib_eigensolver_host ${sources})
target_link_libraries(lib_eigensolver_host pthread)
target_link_libraries(lib_eigensolver_host gmock)
target_link_libraries(lib_eigensolver_host gtest)

add_executable(test_eigensolver_host ../test_main.cpp)
target_link_libraries_whole_archive( test_eigensolver_host lib_eigensolver_host )



set(sources_cuda 
	test_Cardano3x3.cu
	)

cuda_add_library(lib_eigensolver_cuda ${sources_cuda} )
target_link_libraries(lib_eigensolver_cuda pthread)
target_link_libraries(lib_eigensolver_cuda gmock)
target_link_libraries(lib_eigensolver_cuda gtest)

cuda_add_executable(test_eigensolver_cuda ../test_main.cpp)
target_link_libraries_whole_archive( test_eigensolver_cuda lib_eigensolver_cuda )



cuda_add_executable(test_eigensolver ../test_main.cpp)
target_link_libraries_whole_archive( test_eigensolver lib_eigensolver_host )
target_link_libraries_whole_archive( test_eigensolver lib_eigensolver_cuda )
