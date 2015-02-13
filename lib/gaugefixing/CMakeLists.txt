project(gaugefixing)

cmake_minimum_required(VERSION 2.8.3)

set( CUDA_ARCHITECTURE "sm_35" CACHE STRING "CUDA architecture to compile for" )
set( SUN "3" CACHE STRING "Gauge group, 2 = SU(2) or 3 = SU(3) (default)" )
set( FILETYPE "HEADERONLY" CACHE STRING "filetype for gauge fields" )
set( FILETYPEOUT "" CACHE STRING "filetype for gauge fields when writing (leave empty to use same as for reading)" )
option( MIXED_PRECISION "switch on mixed precision (in the SP application)" OFF )
option( USE_TIMESLICEPATTERN "uses a timeslice ordered pattern in LandauGauge" OFF )
option( PRINT_KERNEL_INFO "print kernel infos" OFF )
option( CULGT_ALLOW_C++11 "allows C++11 compilation for cuLGT" ON )

#========================================================


#include( versioning.cmake )

include_directories( $ENV{CULGT_HOME}/include )
link_directories( /usr/local/cuda/lib64)
add_definitions(-std=c++11)
set(CMAKE_CXX_FLAGS "-Wall -O3")

find_package(CUDA)
include(FindCUDA)

if( ${PRINT_KERNEL_INFO} )
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -Xptxas -v" )
endif()

if( ${MIXED_PRECISION} )
    set( CULGT_OPTIONS "${CULGT_OPTIONS} -DDPUPDATE" )
endif()

if( ${USE_TIMESLICEPATTERN} )
    set( CULGT_OPTIONS "${CULGT_OPTIONS} -DCULGT_USE_TIMESLICE_PATTERN" ) 
endif()

set( CULGT_OPTIONS "${CULGT_OPTIONS} -DCULGT_SUN=${SUN}" ) 

if( FILETYPE STREQUAL "VOGT" )
    set( CULGT_OPTIONS "${CULGT_OPTIONS} -DCULGT_FILETYPE_VOGT" ) 
elseif( FILETYPE STREQUAL "HIREP" )
    set( CULGT_OPTIONS "${CULGT_OPTIONS} -DCULGT_FILETYPE_HIREP" )
elseif( FILETYPE STREQUAL "HEADERONLY" )
else()
     message( FATAL_ERROR "FILETYPE " ${FILETYPE} " unknown!" )
endif()

if( FILETYPEOUT STREQUAL "VOGT" )
    set( CULGT_OPTIONS "${CULGT_OPTIONS} -DCULGT_FILETYPE_VOGT_OUT" ) 
elseif( FILETYPEOUT STREQUAL "HIREP" )
    set( CULGT_OPTIONS "${CULGT_OPTIONS} -DCULGT_FILETYPE_HIREP_OUT" )
elseif( FILETYPEOUT STREQUAL "HEADERONLY" )
    set( CULGT_OPTIONS "${CULGT_OPTIONS} -DCULGT_FILETYPE_HEADERONLY_OUT" )
elseif( FILETYPEOUT STREQUAL "" )
else()
     message( FATAL_ERROR "FILETYPEOUT " ${FILETYPEOUT} " unknown!" )
endif()

if( ${CUDA_VERSION} VERSION_LESS 6.5 OR NOT CULGT_ALLOW_C++11 )
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=${CUDA_ARCHITECTURE} -Xptxas -dlcm=cg --use_fast_math")
    message(STATUS "Cannot use -std=c++11 (CUDA version ${CUDA_VERSION})")
else()
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11 -arch=${CUDA_ARCHITECTURE} -Xptxas -dlcm=cg --use_fast_math")
    message(STATUS "Using -std=c++11 (CUDA version ${CUDA_VERSION})")
endif()

#set( CMAKE_VERBOSE_MAKEFILE on )

add_custom_target( MakeGitVersion 
    COMMAND ${CMAKE_COMMAND} -DPATH_TO_CURRENT_PROJECT=${CMAKE_CURRENT_LIST_DIR} -P $ENV{CULGT_HOME}/cmake/versioning.cmake )

cuda_add_executable(CoulombGaugeSP CoulombGaugeFixingApp.cu $ENV{CULGT_HOME}/include/cudacommon/DeviceProperties.cu $ENV{CULGT_HOME}/include/cuLGT1legacy/Chronotimer.cc OPTIONS ${CULGT_OPTIONS})
add_dependencies( CoulombGaugeSP MakeGitVersion )
target_link_libraries(CoulombGaugeSP pthread)
target_link_libraries(CoulombGaugeSP cudart)
target_link_libraries(CoulombGaugeSP cufft)
target_link_libraries(CoulombGaugeSP boost_program_options)

cuda_add_executable(CoulombGaugeDP CoulombGaugeFixingApp.cu $ENV{CULGT_HOME}/include/cudacommon/DeviceProperties.cu $ENV{CULGT_HOME}/include/cuLGT1legacy/Chronotimer.cc OPTIONS -DDOUBLEPRECISION ${CULGT_OPTIONS})
add_dependencies( CoulombGaugeDP MakeGitVersion )
target_link_libraries(CoulombGaugeDP pthread)
target_link_libraries(CoulombGaugeDP cudart)
target_link_libraries(CoulombGaugeDP cufft)
target_link_libraries(CoulombGaugeDP boost_program_options)

cuda_add_executable(LandauGaugeSP LandauGaugeFixingApp.cu $ENV{CULGT_HOME}/include/cudacommon/DeviceProperties.cu $ENV{CULGT_HOME}/include/cuLGT1legacy/Chronotimer.cc OPTIONS ${CULGT_OPTIONS} )
add_dependencies( LandauGaugeSP MakeGitVersion )
target_link_libraries(LandauGaugeSP pthread)
target_link_libraries(LandauGaugeSP cudart)
target_link_libraries(LandauGaugeSP cufft)
target_link_libraries(LandauGaugeSP boost_program_options)

cuda_add_executable(LandauGaugeDP LandauGaugeFixingApp.cu $ENV{CULGT_HOME}/include/cudacommon/DeviceProperties.cu $ENV{CULGT_HOME}/include/cuLGT1legacy/Chronotimer.cc OPTIONS -DDOUBLEPRECISION ${CULGT_OPTIONS})
add_dependencies( LandauGaugeDP MakeGitVersion )
target_link_libraries(LandauGaugeDP pthread)
target_link_libraries(LandauGaugeDP cudart)
target_link_libraries(LandauGaugeDP cufft)
target_link_libraries(LandauGaugeDP boost_program_options)


