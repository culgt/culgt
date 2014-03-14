project(test_lattice)

cmake_minimum_required(VERSION 2.8)

include( ../commonsettings.cmake )

#set( CMAKE_VERBOSE_MAKEFILE on )

set(sources 
	test_Plaquette.cc
	)
add_executable(test_observables_host ../test_main.cpp ${sources})
target_link_libraries(test_observables_host pthread)
target_link_libraries(test_observables_host gmock)
target_link_libraries(test_observables_host gtest)

if( DEFINED allsources )
addWithDirectory( observables "${sources}" )
endif()



set(sources_cuda 
	test_PlaquetteAverage.cu
	)

cuda_add_executable(test_observables_cuda ../test_main.cpp ${sources_cuda} )
target_link_libraries(test_observables_cuda pthread)
target_link_libraries(test_observables_cuda gmock)
target_link_libraries(test_observables_cuda gtest)
target_link_libraries(test_observables_cuda cudart)

if( DEFINED allsources )
addWithDirectory( observables "${sources_cuda}"  )
endif()



cuda_add_executable(test_observables ../test_main.cpp ${sources} ${sources_cuda})
target_link_libraries(test_observables pthread)
target_link_libraries(test_observables gmock)
target_link_libraries(test_observables gtest)



