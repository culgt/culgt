project( culgt_cudacommon )

cmake_minimum_required( VERSION 2.8 )
	
find_package(CUDA)
include(FindCUDA)

set(CMAKE_CXX_FLAGS "-Wall")

cuda_add_library( culgt_cudacommon 
    ../../include/cudacommon/DeviceMemoryManager.cu
    ../../include/cudacommon/DeviceProperties.cu
    )