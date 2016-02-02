project(cuLGT4MILC)

cmake_minimum_required(VERSION 2.8.3)

set( CUDA_ARCHITECTURE "sm_35" CACHE STRING "CUDA architecture to compile for" )
set( CULGT_HOME $ENV{CULGT_HOME} CACHE FILEPATH "path to cuLGT root directory" )
option( DISABLE_AUTOTUNE "disables the autotuner (for debugging purpose)" OFF )
set( CULGT_PATH_TO_MILC $ENV{CULGT_PATH_TO_MILC} CACHE FILEPATH "path to MILC code" )
set( CULGT_PATH_TO_MILC_PROJECT $ENV{CULGT_PATH_TO_MILC_PROJECT} CACHE FILEPATH "path to MILC project directory" )

#========================================================

if( CULGT_HOME STREQUAL "" )
    message( WARNING "Did you set the CULGT_HOME variable?" )
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CULGT_HOME}/cmake/ )

include_directories( ${CULGT_HOME}/include )

include_directories( ${CULGT_HOME}/include )
link_directories( /usr/local/cuda/lib64)
add_definitions(-std=c++11)
set(CMAKE_CXX_FLAGS "-Wall -O3")

if( CULGT_PATH_TO_MILC STREQUAL "" )
    message( FATAL_ERROR "Please set CULGT_PATH_TO_MILC (in cmake or as environmental variable)" )
endif()
if( CULGT_PATH_TO_MILC_PROJECT STREQUAL "" )
    message( FATAL_ERROR "Please set CULGT_PATH_TO_MILC_PROJECT (in cmake or as environmental variable)" )
endif()

include_directories( ${CULGT_PATH_TO_MILC}/include )
include_directories( ${CULGT_PATH_TO_MILC_PROJECT} )

find_package(CUDA)

set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=${CUDA_ARCHITECTURE} -Xptxas -v -Xptxas -dlcm=cg --use_fast_math")

cuda_add_library(MILCcuLGT_SP SHARED MILCcuLGT.cu
    ${CULGT_HOME}/include/cudacommon/DeviceProperties.cu
    ${CULGT_HOME}/include/cudacommon/DeviceMemoryManager.cu
    ${CULGT_HOME}/include/util/timer/Chronotimer.cc)
target_link_libraries(MILCcuLGT_SP pthread)
