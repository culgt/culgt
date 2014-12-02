project(test_lattice)

cmake_minimum_required(VERSION 2.6)

include_directories( ../../include )
include( ../commonsettings.cmake )


#set( CMAKE_VERBOSE_MAKEFILE on )

set(sources 
	test_PowerIterator.cc
	test_Cardano3x3.cc
	test_Matrix.cc
	)
add_executable(test_eigensolver_host ../test_main.cpp ${sources})
target_link_libraries(test_eigensolver_host pthread)
target_link_libraries(test_eigensolver_host gmock)
target_link_libraries(test_eigensolver_host gtest)

if( DEFINED allsources )
addWithDirectory( eigensolver "${sources}" )
endif()



set(sources_cuda 
	test_Cardano3x3.cu
	)

cuda_add_executable(test_eigensolver_cuda ../test_main.cpp ${sources_cuda} )
target_link_libraries(test_eigensolver_cuda pthread)
target_link_libraries(test_eigensolver_cuda gmock)
target_link_libraries(test_eigensolver_cuda gtest)
target_link_libraries(test_eigensolver_cuda cudart)

if( DEFINED allsources )
addWithDirectory( eigensolver "${sources_cuda}"  )
endif()



cuda_add_executable(test_eigensolver ../test_main.cpp ${sources} ${sources_cuda})
target_link_libraries(test_eigensolver pthread)
target_link_libraries(test_eigensolver gmock)
target_link_libraries(test_eigensolver gtest)





