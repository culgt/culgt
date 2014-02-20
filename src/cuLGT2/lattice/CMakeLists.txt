project(GaugeConfiguration)
cmake_minimum_required(VERSION 2.8)

include_directories($ENV{GMOCK_HOME}/include $ENV{GMOCK_HOME}/gtest/include)
link_directories($ENV{GMOCK_HOME}/mybuild $ENV{GMOCK_HOME}/gtest/mybuild)
add_definitions()
set(CMAKE_CXX_FLAGS "${CMAXE_CXX_FLAGS} -Wall")

set(sources 
	main.cpp 
	LocalLink_test.cc
	su3/ParameterizationMediatorSU3_Real12_Real18_test.cc 
	su3/SU3.cc
	configuration_patterns/test_StandardPattern.cc
	configuration_patterns/GPUPattern.cc
	configuration_patterns/GPUPatternParityPriority.cc
	GlobalLink.cc 
	)
add_executable(test ${sources})
target_link_libraries(test pthread)
target_link_libraries(test gmock)
target_link_libraries(test gtest)

set( sources_StandardPattern
	main.cpp 
	configuration_patterns/test_StandardPattern.cc
	)
add_executable(test_standard_pattern ${sources_StandardPattern}	)
target_link_libraries(test_standard_pattern pthread)
target_link_libraries(test_standard_pattern gmock)
target_link_libraries(test_standard_pattern gtest)

