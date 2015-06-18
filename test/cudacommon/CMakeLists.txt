project(test_cudacommon)
cmake_minimum_required(VERSION 2.8)

include( ../commonsettings.cmake )


#set(sources 
#	)
#	
#add_library(lib_cudacommon_host ${sources})
#target_link_libraries(lib_cudacommon_host pthread)
#target_link_libraries(lib_cudacommon_host gmock)
#target_link_libraries(lib_cudacommon_host gtest)
#target_link_libraries(lib_cudacommon_host boost_program_options)
#
#add_executable( test_cudacommon_host ../test_main.cpp )
#target_link_libraries_whole_archive( test_cudacommon_host lib_cudacommon_host )



set(sources_cuda 
	test_DeviceCommunicator.cu
	)

cuda_add_library(lib_cudacommon_cuda ${sources_cuda} )
target_link_libraries(lib_cudacommon_cuda pthread)
target_link_libraries(lib_cudacommon_cuda gmock)
target_link_libraries(lib_cudacommon_cuda gtest)
target_link_libraries(lib_cudacommon_cuda boost_program_options)

cuda_add_executable( test_cudacommon_cuda ../test_main.cpp )
target_link_libraries_whole_archive( test_cudacommon_cuda lib_cudacommon_cuda )



cuda_add_executable(test_cudacommon ../test_main.cpp )
#target_link_libraries_whole_archive( test_cudacommon lib_cudacommon_host )
target_link_libraries_whole_archive( test_cudacommon lib_cudacommon_cuda )
