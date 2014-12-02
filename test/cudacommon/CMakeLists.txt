project(test_cudacommon)
cmake_minimum_required(VERSION 2.8)

include( ../commonsettings.cmake )
include_directories(../../include)

set(sources 
	
	)
add_executable(test_cudacommon_host ../test_main.cpp  ${sources})
target_link_libraries(test_cudacommon_host pthread)
target_link_libraries(test_cudacommon_host gmock)
target_link_libraries(test_cudacommon_host gtest)
target_link_libraries(test_cudacommon_host boost_program_options)

if( DEFINED allsources )
addWithDirectory( cudacommon "${sources}" )
endif()

set(sources_cuda 
	test_DeviceCommunicator.cu
	)

cuda_add_executable(test_cudacommon_cuda ../test_main.cpp ${sources_cuda} )
target_link_libraries(test_cudacommon_cuda pthread)
target_link_libraries(test_cudacommon_cuda gmock)
target_link_libraries(test_cudacommon_cuda gtest)
target_link_libraries(test_cudacommon_cuda cudart)
target_link_libraries(test_cudacommon_cuda boost_program_options)

if( DEFINED allsources )
addWithDirectory( cudacommon "${sources_cuda}"  )
endif()




cuda_add_executable(test_cudacommon ../test_main.cpp ${sources} ${sources_cuda})
target_link_libraries(test_cudacommon pthread)
target_link_libraries(test_cudacommon gmock)
target_link_libraries(test_cudacommon gtest)
target_link_libraries(test_cudacommon boost_program_options)