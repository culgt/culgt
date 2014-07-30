project(test_application)
cmake_minimum_required(VERSION 2.8)

include( ../commonsettings.cmake )

set(sources 
	test_FileIterator.cc
	)
add_executable(test_application_host ../test_main.cpp  ${sources})
target_link_libraries(test_application_host pthread)
target_link_libraries(test_application_host gmock)
target_link_libraries(test_application_host gtest)
target_link_libraries(test_application_host boost_program_options)

if( DEFINED allsources )
addWithDirectory( application "${sources}" )
endif()

set(sources_cuda 
	test_GaugeConfigurationIteratingApplication.cc
	)

cuda_add_executable(test_application_cuda ../test_main.cpp ${sources_cuda} )
target_link_libraries(test_application_cuda pthread)
target_link_libraries(test_application_cuda gmock)
target_link_libraries(test_application_cuda gtest)
target_link_libraries(test_application_cuda cudart)
target_link_libraries(test_application_cuda boost_program_options)

if( DEFINED allsources )
addWithDirectory( application "${sources_cuda}"  )
endif()




cuda_add_executable(test_application ../test_main.cpp ${sources} ${sources_cuda})
target_link_libraries(test_application pthread)
target_link_libraries(test_application gmock)
target_link_libraries(test_application gtest)
target_link_libraries(test_application boost_program_options)