project(test_math)
cmake_minimum_required(VERSION 2.8)

include_directories(../../include)
include( ../commonsettings.cmake )

set(sources 
	)
add_executable(test_math_host ../test_main.cpp  ${sources})
target_link_libraries(test_math_host pthread)
target_link_libraries(test_math_host gmock)
target_link_libraries(test_math_host gtest)

if( DEFINED allsources )
addWithDirectory( math "${sources}" )
endif()

set(sources_cuda 
	test_Complex.cu
	)

cuda_add_executable(test_math_cuda ../test_main.cpp ${sources_cuda} )
target_link_libraries(test_math_cuda pthread)
target_link_libraries(test_math_cuda gmock)
target_link_libraries(test_math_cuda gtest)
target_link_libraries(test_math_cuda cudart)

if( DEFINED allsources )
addWithDirectory( math "${sources_cuda}"  )
endif()




cuda_add_executable(test_math ../test_main.cpp ${sources} ${sources_cuda})
target_link_libraries(test_math pthread)
target_link_libraries(test_math gmock)
target_link_libraries(test_math gtest)