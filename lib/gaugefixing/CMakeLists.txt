project(propagator)

cmake_minimum_required(VERSION 2.8)


include_directories( $ENV{CULGT_HOME}/include )
include_directories( ../../include )
link_directories( /usr/local/cuda/lib64)
add_definitions(-std=c++11)
set(CMAKE_CXX_FLAGS "-Wall -O3")

find_package(CUDA)
include(FindCUDA)

if( ${CUDA_VERSION} VERSION_LESS 6.5 )
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=sm_20 -Xptxas -v -Xptxas -dlcm=cg --use_fast_math")
    message(STATUS "Cannot use -std=c++11 (CUDA version ${CUDA_VERSION})")
else()
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11 -arch=sm_20 -Xptxas -v -Xptxas -dlcm=cg --use_fast_math")
    message(STATUS "Using -std=c++11 (CUDA version ${CUDA_VERSION})")
endif()

#set( CMAKE_VERBOSE_MAKEFILE on )

cuda_add_executable(CoulombGaugeDP CoulombGaugeFixingApp.cu $ENV{CULGT_HOME}/include/cuLGT1legacy/Chronotimer.cc OPTIONS -DDOUBLEPRECISION)
target_link_libraries(CoulombGaugeDP pthread)
target_link_libraries(CoulombGaugeDP cudart)
target_link_libraries(CoulombGaugeDP cufft)
target_link_libraries(CoulombGaugeDP boost_program_options)

cuda_add_executable(CoulombGaugeSP CoulombGaugeFixingApp.cu $ENV{CULGT_HOME}/include/cuLGT1legacy/Chronotimer.cc)
target_link_libraries(CoulombGaugeSP pthread)
target_link_libraries(CoulombGaugeSP cudart)
target_link_libraries(CoulombGaugeSP cufft)
target_link_libraries(CoulombGaugeSP boost_program_options)

cuda_add_executable(LandauGaugeSP LandauGaugeFixingApp.cu $ENV{CULGT_HOME}/include/cuLGT1legacy/Chronotimer.cc)
target_link_libraries(LandauGaugeSP pthread)
target_link_libraries(LandauGaugeSP cudart)
target_link_libraries(LandauGaugeSP cufft)
target_link_libraries(LandauGaugeSP boost_program_options)

cuda_add_executable(LandauGaugeDP LandauGaugeFixingApp.cu $ENV{CULGT_HOME}/include/cuLGT1legacy/Chronotimer.cc OPTIONS -DDOUBLEPRECISION)
target_link_libraries(LandauGaugeDP pthread)
target_link_libraries(LandauGaugeDP cudart)
target_link_libraries(LandauGaugeDP cufft)
target_link_libraries(LandauGaugeDP boost_program_options)


