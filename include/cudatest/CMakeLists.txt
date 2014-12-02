project(cuda_test_helper_test)

#The FindCUDA script is distributed since version 2.8
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
 
add_definitions()
set(CMAKE_CXX_FLAGS "${CMAXE_CXX_FLAGS} -Wall")
include_directories($ENV{GMOCK_HOME}/include $ENV{GMOCK_HOME}/gtest/include)
link_directories($ENV{GMOCK_HOME}/mybuild $ENV{GMOCK_HOME}/gtest/mybuild)
 
 
 set(sources_host 
	main.cpp
	test_example_gtest_plugin_host.cc
	)
add_executable(test_host ${sources_host} )
target_link_libraries(test_host pthread)
target_link_libraries(test_host gmock)
target_link_libraries(test_host gtest)
 
FIND_PACKAGE(CUDA)
 
INCLUDE(FindCUDA)


set(sources 
	main.cpp
	test_example_gtest_plugin.cu
	)



set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-arch=sm_20;--ptxas-options=-v;)

CUDA_ADD_EXECUTABLE(test ${sources} )
target_link_libraries(test pthread)
target_link_libraries(test gmock)
target_link_libraries(test gtest)


	

