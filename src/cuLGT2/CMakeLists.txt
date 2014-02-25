project(test_main)
cmake_minimum_required(VERSION 2.8)

include( addWithDirectory.cmake )

include_directories($ENV{GMOCK_HOME}/include $ENV{GMOCK_HOME}/gtest/include)
link_directories($ENV{GMOCK_HOME}/mybuild $ENV{GMOCK_HOME}/gtest/mybuild)
add_definitions()
set(CMAKE_CXX_FLAGS "${CMAXE_CXX_FLAGS} -Wall")

set( allsources "" )

add_subdirectory( lattice )
add_subdirectory( cuLGT1legacy )
add_subdirectory( propagators )

find_package(CUDA)
include(FindCUDA)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-arch=sm_20;)

message( "allsources = ${allsources}" )

CUDA_ADD_EXECUTABLE(test test_main.cpp ${allsources} )
target_link_libraries(test pthread)
target_link_libraries(test gmock)
target_link_libraries(test gtest)
