project(test_main)
cmake_minimum_required(VERSION 2.8)

include( commonsettings.cmake )
include( addWithDirectory.cmake )

include_directories( ../include )

set( allsources "" )

add_subdirectory( application )
add_subdirectory( cudacommon )
add_subdirectory( cuLGT1legacy )
add_subdirectory( gaugefixing )
add_subdirectory( lattice )
add_subdirectory( eigensolver )

message( "allsources = ${allsources}" )

CUDA_ADD_EXECUTABLE(test test_main.cpp ${allsources} $ENV{CULGT_HOME}/include/cudacommon/DeviceProperties.cu )
target_link_libraries(test pthread)
target_link_libraries(test gmock)
target_link_libraries(test gtest)
target_link_libraries(test cufft)
target_link_libraries(test boost_program_options)
