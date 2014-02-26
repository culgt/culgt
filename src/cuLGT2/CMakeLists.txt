project(test_main)
cmake_minimum_required(VERSION 2.8)

include( commonsettings.cmake )
include( addWithDirectory.cmake )

set( allsources "" )

add_subdirectory( lattice )
add_subdirectory( cuLGT1legacy )
add_subdirectory( propagators )

message( "allsources = ${allsources}" )

CUDA_ADD_EXECUTABLE(test test_main.cpp ${allsources} )
target_link_libraries(test pthread)
target_link_libraries(test gmock)
target_link_libraries(test gtest)
