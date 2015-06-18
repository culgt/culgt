project(test)
cmake_minimum_required(VERSION 2.8)

include( commonsettings.cmake )


add_subdirectory( application )
add_subdirectory( cudacommon )
add_subdirectory( eigensolver )
add_subdirectory( lattice )
add_subdirectory( math )
add_subdirectory( observables )
add_subdirectory( util )


add_executable( test_host test_main.cpp )
target_link_libraries_whole_archive( test_host lib_application_host )
#target_link_libraries_whole_archive( test_host lib_cudacommon_host )
target_link_libraries_whole_archive( test_host lib_eigensolver_host )
target_link_libraries_whole_archive( test_host lib_lattice_host )
#target_link_libraries_whole_archive( test_host lib_math_host )
target_link_libraries_whole_archive( test_host lib_observables_host )
target_link_libraries_whole_archive( test_host lib_util_host )



cuda_add_executable( test_cuda test_main.cpp )
target_link_libraries_whole_archive( test_cuda lib_application_cuda )
target_link_libraries_whole_archive( test_cuda lib_cudacommon_cuda )
target_link_libraries_whole_archive( test_cuda lib_eigensolver_cuda )
target_link_libraries_whole_archive( test_cuda lib_lattice_cuda )
target_link_libraries_whole_archive( test_cuda lib_math_cuda )
target_link_libraries_whole_archive( test_cuda lib_observables_cuda )
target_link_libraries_whole_archive( test_cuda lib_util_cuda )



cuda_add_executable( test_all test_main.cpp )
target_link_libraries_whole_archive( test_all lib_application_host )
target_link_libraries_whole_archive( test_all lib_application_cuda )
target_link_libraries_whole_archive( test_all lib_cudacommon_cuda )
target_link_libraries_whole_archive( test_all lib_eigensolver_host )
target_link_libraries_whole_archive( test_all lib_eigensolver_cuda )
target_link_libraries_whole_archive( test_all lib_lattice_host )
target_link_libraries_whole_archive( test_all lib_lattice_cuda )
target_link_libraries_whole_archive( test_all lib_math_cuda )
target_link_libraries_whole_archive( test_all lib_observables_host )
target_link_libraries_whole_archive( test_all lib_observables_cuda )
target_link_libraries_whole_archive( test_all lib_util_host )
target_link_libraries_whole_archive( test_all lib_util_cuda )

