project(GaugeConfiguration)
cmake_minimum_required(VERSION 2.8)

include_directories($ENV{GMOCK_HOME}/include $ENV{GMOCK_HOME}/gtest/include)
link_directories($ENV{GMOCK_HOME}/mybuild $ENV{GMOCK_HOME}/gtest/mybuild)
add_definitions()
set(CMAKE_CXX_FLAGS "${CMAXE_CXX_FLAGS} -Wall")

set(sources 
	main.cpp 
	test_LocalLink.cc
	su3/ParameterizationMediatorSU3_Real12_Real18_test.cc 
	su3/SU3.cc
	configuration_patterns/test_StandardPattern.cc
	configuration_patterns/GPUPattern.cc
	configuration_patterns/GPUPatternParityPriority.cc
	test_GlobalLink.cc 
	test_GlobalLinkWithPatterns.cc
	test_LocalGlobalLinkCompatibility.cc
	)
add_executable(test ${sources})
target_link_libraries(test pthread)
target_link_libraries(test gmock)
target_link_libraries(test gtest)


set(sources2 
	main.cpp 
	test_LocalGlobalLinkCompatibility.cc
	)
add_executable(test2 ${sources2})
target_link_libraries(test2 pthread)
target_link_libraries(test2 gmock)
target_link_libraries(test2 gtest)
