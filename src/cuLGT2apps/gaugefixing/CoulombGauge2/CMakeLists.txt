project(propagator)

cmake_minimum_required(VERSION 2.8)


include_directories( ../../../cuLGT2 )
link_directories( /usr/local/cuda/lib64)
add_definitions(-std=c++11)
set(CMAKE_CXX_FLAGS "-Wall -O3")

find_package(CUDA)
include(FindCUDA)
set(CUDA_NVCC_FLAGS "-arch=sm_20 -Xptxas -v")

#set( CMAKE_VERBOSE_MAKEFILE on )

cuda_add_executable(DMCGaugeDP DirectMaximalCenterGaugeFixingApp.cu ../../../cuLGT2/cuLGT1legacy/Chronotimer.cc OPTIONS -DDOUBLEPRECISION )
target_link_libraries(DMCGaugeDP pthread)
target_link_libraries(DMCGaugeDP cudart)
target_link_libraries(DMCGaugeDP cufft)
target_link_libraries(DMCGaugeDP boost_program_options)

cuda_add_executable(DMCGaugeSP DirectMaximalCenterGaugeFixingApp.cu ../../../cuLGT2/cuLGT1legacy/Chronotimer.cc)
target_link_libraries(DMCGaugeSP pthread)
target_link_libraries(DMCGaugeSP cudart)
target_link_libraries(DMCGaugeSP cufft)
target_link_libraries(DMCGaugeSP boost_program_options)