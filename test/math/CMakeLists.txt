project(test_math)
cmake_minimum_required(VERSION 2.8)

include( ../commonsettings.cmake )



set(sources_cuda 
	test_Complex.cu
	)

cuda_add_library(lib_math_cuda ${sources_cuda} )
target_link_libraries(lib_math_cuda pthread)
target_link_libraries(lib_math_cuda gmock)
target_link_libraries(lib_math_cuda gtest)

cuda_add_executable( test_math_cuda ../test_main.cpp )
target_link_libraries_whole_archive( test_math_cuda lib_math_cuda )



cuda_add_executable(test_math ../test_main.cpp )
target_link_libraries_whole_archive( test_math lib_math_cuda )
