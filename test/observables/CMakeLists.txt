project(test_observables)
cmake_minimum_required(VERSION 2.8)

include( ../commonsettings.cmake )



set(sources 
	test_Plaquette.cc
	)
add_library(lib_observables_host  ${sources})
target_link_libraries(lib_observables_host pthread)
target_link_libraries(lib_observables_host gmock)
target_link_libraries(lib_observables_host gtest)

add_executable( test_observables_host ../test_main.cpp )
target_link_libraries_whole_archive( test_observables_host lib_observables_host )



set(sources_cuda
	test_PlaquetteAverage.cu
	)

cuda_add_library(lib_observables_cuda ${sources_cuda} )
target_link_libraries(lib_observables_cuda pthread)
target_link_libraries(lib_observables_cuda gmock)
target_link_libraries(lib_observables_cuda gtest)
target_link_libraries(lib_observables_cuda culgt_cudacommon )

cuda_add_executable(test_observables_cuda ../test_main.cpp )
target_link_libraries_whole_archive( test_observables_cuda lib_observables_cuda )



cuda_add_executable(test_observables ../test_main.cpp )
target_link_libraries_whole_archive( test_observables lib_observables_host )
target_link_libraries_whole_archive( test_observables lib_observables_cuda )



