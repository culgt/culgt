project(test_application)
cmake_minimum_required(VERSION 2.8)

include( ../commonsettings.cmake )



set(sources 
	test_FileIterator.cc
	)
	
add_library( lib_application_host ${sources} )
target_link_libraries( lib_application_host pthread)
target_link_libraries( lib_application_host gmock)
target_link_libraries( lib_application_host gtest)
target_link_libraries( lib_application_host liblime.a)
	
add_executable( test_application_host ../test_main.cpp )
target_link_libraries_whole_archive( test_application_host lib_application_host )



set(sources_cuda 
	test_GaugeConfigurationIteratingApplication.cu
	)

cuda_add_library(lib_application_cuda ${sources_cuda} )
target_link_libraries(lib_application_cuda pthread)
target_link_libraries(lib_application_cuda gmock)
target_link_libraries(lib_application_cuda gtest)
target_link_libraries(lib_application_cuda boost_program_options)
target_link_libraries(lib_application_cuda culgt_cudacommon )
target_link_libraries( lib_application_cuda liblime.a)
target_link_libraries( lib_application_cuda tinyxml)

cuda_add_executable( test_application_cuda ../test_main.cpp )
target_link_libraries_whole_archive( test_application_cuda lib_application_cuda )



cuda_add_executable(test_application ../test_main.cpp )
target_link_libraries_whole_archive( test_application lib_application_host )
target_link_libraries_whole_archive( test_application lib_application_cuda )
