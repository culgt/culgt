function( makeVersionString myString )
    set( NEW_VERSION_STRING  "#define CULGT_VERSION_MACRO \"${myString}\"" PARENT_SCOPE )
endfunction( makeVersionString myString ) 

function( writeVersionString )
    execute_process( COMMAND bash -c "echo '${NEW_VERSION_STRING}' > ${PATH_TO_CURRENT_PROJECT}/version" )
endfunction( writeVersionString )

# check if we are in a git managed project
execute_process( COMMAND git branch
                RESULT_VARIABLE GIT_EXIT_CODE 
                OUTPUT_QUIET ERROR_QUIET )
if( GIT_EXIT_CODE EQUAL 0 )
    set( IS_GIT_DIRECTORY ON )
else()
    set( IS_GIT_DIRECTORY OFF )
endif()
# now IS_GIT_DIRECTORY is ON or OFF 

if( IS_GIT_DIRECTORY )
    # get version from git
    execute_process(COMMAND git describe --long --dirty
                    OUTPUT_VARIABLE CULGT_VERSION
                    ERROR_VARIABLE DONTPRINTERRORS
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    RESULT_VARIABLE GIT_EXIT_CODE )
    
    if( GIT_EXIT_CODE EQUAL 0 )
        # read old version from version file
        if( EXISTS ${CMAKE_CURRENT_LIST_DIR}/version )
            execute_process( COMMAND cat "${PATH_TO_CURRENT_PROJECT}/version" OUTPUT_VARIABLE OLD_VERSION_STRING OUTPUT_STRIP_TRAILING_WHITESPACE )
        endif( EXISTS ${CMAKE_CURRENT_LIST_DIR}/version )              
                        
        makeVersionString( ${CULGT_VERSION} )                
        
        # write new version file only if different from old version (avoids recompile if nothing changed)
        if( NOT OLD_VERSION_STRING STREQUAL NEW_VERSION_STRING )
            writeVersionString()
        endif()
    
    else()
        message( "No cuLGT version number available" )
        makeVersionString( "undefined" )       
        writeVersionString() 
    endif()
    
else( IS_GIT_DIRECTORY )
    # if we are not git managed, a version file should come with the distribution (if not set version to undefined)
    if( NOT EXISTS ${CMAKE_CURRENT_LIST_DIR}/version )
        message( "No cuLGT version number available" )
        makeVersionString( "undefined" )       
        writeVersionString() 
    endif( NOT EXISTS ${CMAKE_CURRENT_LIST_DIR}/version )
endif( IS_GIT_DIRECTORY )