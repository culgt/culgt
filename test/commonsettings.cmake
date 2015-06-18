if( NOT DEFINED commonsettings )
set( commonsettings "")

    set( CULGT_HOME $ENV{CULGT_HOME} )
    if( NOT DEFINED CULGT_HOME )
        message( FATAL_ERROR "Please set CULGT_HOME (in cmake or as environmental variable)" )
    endif()
    
    include( ${CULGT_HOME}/test/gmock_setup.cmake )
    include( ${CULGT_HOME}/cmake/target_link_libraries_whole_archive.cmake )

    set( CUDA_ARCH "sm_35" CACHE STRING "CUDA architecture to compile for" )
    option( PRINT_KERNEL_INFO "print kernel infos" OFF )
    option( CULGT_ALLOW_C++11 "allows C++11 compilation for cuLGT" ON )
    



#    set( CMAKE_VERBOSE_MAKEFILE on )
    
	include_directories(
	    ${GMOCK_HOME}/include
	    ${GTEST_HOME}/include
	    ${CULGT_HOME}/include
	    )
	link_directories(
	    ${GMOCK_LIBDIR}
	    ${GTEST_LIBDIR}
	    )
	    
	add_definitions(-std=c++11)
	set(CMAKE_CXX_FLAGS "-Wall")
	
	find_package(CUDA)
	include(FindCUDA)
	
    message(STATUS "Compiling for cuda architecture ${CUDA_ARCH}" )
	set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-arch=${CUDA_ARCH}" )
	
	if( ${PRINT_KERNEL_INFO} )
        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -Xptxas -v" )
    endif()
    
    
    if( ${CUDA_VERSION} VERSION_LESS 6.5 )
        message( STATUS "Cannot use c++11 in nvcc code (cuda version <= 6.5)" )
    elseif( ${CULGT_ALLOW_C++11} )
	    message( STATUS "Using c++11 in nvcc code" )
        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11" )
    else()
        message( STATUS "c++11 in nvcc code disabled by user" )
    endif()
	

    add_subdirectory( ${CULGT_HOME}/lib/cudacommon culgt_cudacommon )

#   copy/create testdata
    file( COPY ${CULGT_HOME}/test/testdata/test_configSU2N4T8SP.vogt DESTINATION "." )
    file( WRITE ${CMAKE_BINARY_DIR}/tune_optimalid "0 0" )

endif()
