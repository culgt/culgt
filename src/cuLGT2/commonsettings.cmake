if( NOT DEFINED commonsettings )
	include_directories($ENV{GMOCK_HOME}/include $ENV{GMOCK_HOME}/gtest/include)
	link_directories($ENV{GMOCK_HOME}/mybuild $ENV{GMOCK_HOME}/gtest/mybuild /usr/local/cuda/lib64)
	add_definitions(-std=c++11)
	set(CMAKE_CXX_FLAGS "-Wall")
	
	find_package(CUDA)
	include(FindCUDA)
	set(CUDA_NVCC_FLAGS "-arch=sm_20 -Xptxas -v")
#    -std=c++11 #available starting with cuda 6.5
endif()
