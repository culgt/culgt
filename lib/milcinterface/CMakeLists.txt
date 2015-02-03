project(propagator)

cmake_minimum_required(VERSION 2.8)


include_directories( $ENV{CULGT_HOME}/include )
include_directories( $ENV{CULGT_PATH_TO_MILC}/include )
include_directories( $ENV{CULGT_PATH_TO_MILC_PROJECT} )
link_directories( /usr/local/cuda/lib64)
set(CMAKE_CXX_FLAGS "-Wall -O3")

find_package(CUDA)
include(FindCUDA)

#if( ${CUDA_VERSION} VERSION_LESS 6.5 )
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=sm_20 -Xptxas -v -Xptxas -dlcm=cg --use_fast_math")
    message(STATUS "Cannot use -std=c++11 (CUDA version ${CUDA_VERSION})")
#else()
#    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11 -arch=sm_20 -Xptxas -v -Xptxas -dlcm=cg --use_fast_math")
#    message(STATUS "Using -std=c++11 (CUDA version ${CUDA_VERSION})")
#endif()

#set( CMAKE_VERBOSE_MAKEFILE on )

cuda_add_library(MILCcuLGT_SP SHARED MILCcuLGT.cu $ENV{CULGT_HOME}/include/cudacommon/DeviceProperties.cu $ENV{CULGT_HOME}/include/cuLGT1legacy/Chronotimer.cc )
target_link_libraries(MILCcuLGT_SP pthread)
target_link_libraries(MILCcuLGT_SP cudart)
target_link_libraries(MILCcuLGT_SP cufft)
target_link_libraries(MILCcuLGT_SP boost_program_options)