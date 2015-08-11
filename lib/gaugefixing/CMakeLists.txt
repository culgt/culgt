project(gaugefixing)

cmake_minimum_required(VERSION 2.8.3)

set( CUDA_ARCHITECTURE "sm_35" CACHE STRING "CUDA architecture to compile for" )
set( SUN "3" CACHE STRING "Gauge group, 2 = SU(2) or 3 = SU(3) (default)" )
option( MIXED_PRECISION "switch on mixed precision (in the SP application)" OFF )
option( USE_TIMESLICEPATTERN "uses a timeslice ordered pattern in LandauGauge" OFF )
option( PRINT_KERNEL_INFO "print kernel infos" OFF )
option( CULGT_ALLOW_C++11 "allows C++11 compilation for cuLGT" ON )
option( DISABLE_AUTOTUNE "disables the autotuner (for debugging purpose)" OFF )

#========================================================

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{CULGT_HOME}/cmake/ )

include_directories( $ENV{CULGT_HOME}/include )
link_directories( /usr/local/cuda/lib64)
add_definitions(-std=c++11)
set(CMAKE_CXX_FLAGS "-Wall -O3")

find_package(CUDA)
include( nvcc_version )

find_package( Boost COMPONENTS program_options REQUIRED )
include_directories( ${Boost_INCLUDE_DIRS} )

include( CheckHaveILDG )

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

if( ${DISABLE_AUTOTUNE} )
    set( CULGT_OPTIONS "${CULGT_OPTIONS} -DCULGT_NO_AUTOTUNE" ) 
endif()

if( ${CUDA_VERSION} VERSION_LESS 6.5 OR NOT CULGT_ALLOW_C++11 )
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=${CUDA_ARCHITECTURE} -Xptxas -dlcm=cg --use_fast_math")
    message(STATUS "Cannot use -std=c++11 (CUDA version ${CUDA_VERSION})")
else()
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11 -arch=${CUDA_ARCHITECTURE} -Xptxas -dlcm=cg --use_fast_math")
    message(STATUS "Using -std=c++11 (CUDA version ${CUDA_VERSION})")
endif()

add_custom_target( MakeGitVersion 
    COMMAND ${CMAKE_COMMAND} -DPATH_TO_CURRENT_PROJECT=${CMAKE_CURRENT_LIST_DIR} -P $ENV{CULGT_HOME}/cmake/versioning.cmake )

cuda_add_library( GaugeFixingLibrary 
    $ENV{CULGT_HOME}/include/cudacommon/DeviceProperties.cu
    $ENV{CULGT_HOME}/include/cudacommon/DeviceMemoryManager.cu
    $ENV{CULGT_HOME}/include/util/timer/Chronotimer.cc
    )
target_link_libraries(GaugeFixingLibrary pthread)
target_link_libraries(GaugeFixingLibrary ${Boost_LIBRARIES} )
target_link_libraries(GaugeFixingLibrary ${CUDA_LIBRARIES} )
target_link_libraries(GaugeFixingLibrary ${lime_LIBRARIES} )
target_link_libraries(GaugeFixingLibrary ${tinyxml_LIBRARIES})

cuda_add_executable(CoulombGaugeSP
    CoulombGaugeFixingApp.cu
    OPTIONS ${CULGT_OPTIONS})
add_dependencies( CoulombGaugeSP MakeGitVersion )
target_link_libraries(CoulombGaugeSP GaugeFixingLibrary )

cuda_add_executable(CoulombGaugeDP
    CoulombGaugeFixingApp.cu
    OPTIONS -DDOUBLEPRECISION ${CULGT_OPTIONS})
add_dependencies( CoulombGaugeDP MakeGitVersion )
target_link_libraries(CoulombGaugeDP GaugeFixingLibrary)

cuda_add_executable(LandauGaugeSP
    LandauGaugeFixingApp.cu
    OPTIONS ${CULGT_OPTIONS} )
add_dependencies( LandauGaugeSP MakeGitVersion )
target_link_libraries(LandauGaugeSP GaugeFixingLibrary)

cuda_add_executable(LandauGaugeDP
    LandauGaugeFixingApp.cu
    OPTIONS -DDOUBLEPRECISION ${CULGT_OPTIONS})
add_dependencies( LandauGaugeDP MakeGitVersion )
target_link_libraries(LandauGaugeDP GaugeFixingLibrary)


