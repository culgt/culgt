
execute_process( COMMAND git branch
                RESULT_VARIABLE GIT_EXIT_CODE 
                OUTPUT_QUIET ERROR_QUIET )
if( GIT_EXIT_CODE EQUAL 0 )
    set( IS_GIT_DIRECTORY ON )
endif()

if( IS_GIT_DIRECTORY )
execute_process(COMMAND git describe --long --dirty
                OUTPUT_VARIABLE CULGT_VERSION
                OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process( COMMAND bash -c "echo \"\#define CULGT_VERSION_MACRO \\\"${CULGT_VERSION}\\\"\" > ${CMAKE_CURRENT_LIST_DIR}/version" )
else( IS_GIT_DIRECTORY )
    if( NOT EXISTS ${CMAKE_SOURCE_DIR}/version )
        message( WARNING "No version number available" )
        execute_process( COMMAND bash -c "echo \"\#define CULGT_VERSION_MACRO \\\"undefined\\\"\" > ${CMAKE_CURRENT_LIST_DIR}/version" )
    endif( NOT EXISTS ${CMAKE_SOURCE_DIR}/version )
endif( IS_GIT_DIRECTORY )