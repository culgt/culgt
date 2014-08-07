project(propagator)

cmake_minimum_required(VERSION 2.8)


include_directories( ../../../cuLGT2 )
link_directories( /usr/local/cuda/lib64)
add_definitions(-std=c++11)
set(CMAKE_CXX_FLAGS "-Wall -O3")

find_package(CUDA)
include(FindCUDA)
set(CUDA_NVCC_FLAGS "-arch=sm_20 -Xptxas -v -Xptxas -dlcm=cg --use_fast_math")

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

cuda_add_executable(CoulombGaugeDP CoulombGaugeFixingApp.cu ../../../cuLGT2/cuLGT1legacy/Chronotimer.cc OPTIONS -DDOUBLEPRECISION)
target_link_libraries(CoulombGaugeDP pthread)
target_link_libraries(CoulombGaugeDP cudart)
target_link_libraries(CoulombGaugeDP cufft)
target_link_libraries(CoulombGaugeDP boost_program_options)

cuda_add_executable(CoulombGaugeSP CoulombGaugeFixingApp.cu ../../../cuLGT2/cuLGT1legacy/Chronotimer.cc)
target_link_libraries(CoulombGaugeSP pthread)
target_link_libraries(CoulombGaugeSP cudart)
target_link_libraries(CoulombGaugeSP cufft)
target_link_libraries(CoulombGaugeSP boost_program_options)

cuda_add_executable(DMCGaugeSpatialDP DirectMaximalCenterGaugeFixingSpatialApp.cu ../../../cuLGT2/cuLGT1legacy/Chronotimer.cc OPTIONS -DDOUBLEPRECISION )
target_link_libraries(DMCGaugeSpatialDP pthread)
target_link_libraries(DMCGaugeSpatialDP cudart)
target_link_libraries(DMCGaugeSpatialDP cufft)
target_link_libraries(DMCGaugeSpatialDP boost_program_options)

cuda_add_executable(DMCGaugeSpatialSP DirectMaximalCenterGaugeFixingSpatialApp.cu ../../../cuLGT2/cuLGT1legacy/Chronotimer.cc)
target_link_libraries(DMCGaugeSpatialSP pthread)
target_link_libraries(DMCGaugeSpatialSP cudart)
target_link_libraries(DMCGaugeSpatialSP cufft)
target_link_libraries(DMCGaugeSpatialSP boost_program_options)