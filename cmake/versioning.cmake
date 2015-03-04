function( makeVersionString myString myVariableName )
    set( NEW_VERSION_STRING  "#define ${myVariableName} \"${myString}\"" PARENT_SCOPE )
endfunction( makeVersionString myString ) 

function( writeVersionString )
    execute_process( COMMAND bash -c "echo '${NEW_VERSION_STRING}' > ${PATH_TO_CURRENT_PROJECT}/${CULGT_VERSION_FILE}" )
endfunction( writeVersionString )

if( NOT DEFINED PATH_TO_FOREIGN_GIT_PROJECT )
    set( PATH_TO_FOREIGN_GIT_PROJECT ${PATH_TO_CURRENT_PROJECT} )
endif()

if( NOT DEFINED CULGT_VERSION_VARIABLE )
    set( CULGT_VERSION_VARIABLE "CULGT_VERSION_MACRO" )
endif()

if( NOT DEFINED CULGT_VERSION_FILE )
    set( CULGT_VERSION_FILE "version" )
endif()

# check if we are in a git managed project
execute_process( COMMAND git branch
                WORKING_DIRECTORY ${PATH_TO_FOREIGN_GIT_PROJECT}
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
                    WORKING_DIRECTORY ${PATH_TO_FOREIGN_GIT_PROJECT}
                    OUTPUT_VARIABLE CULGT_VERSION
                    ERROR_VARIABLE DONTPRINTERRORS
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    RESULT_VARIABLE GIT_EXIT_CODE )
    
    if( GIT_EXIT_CODE EQUAL 0 )
        # read old version from version file
        if( EXISTS ${PATH_TO_CURRENT_PROJECT}/${CULGT_VERSION_FILE} )
            execute_process( COMMAND cat "${PATH_TO_CURRENT_PROJECT}/${CULGT_VERSION_FILE}"
                                WORKING_DIRECTORY ${PATH_TO_FOREIGN_GIT_PROJECT}
                                OUTPUT_VARIABLE OLD_VERSION_STRING OUTPUT_STRIP_TRAILING_WHITESPACE )
        endif( EXISTS ${PATH_TO_CURRENT_PROJECT}/${CULGT_VERSION_FILE} )              
                       
        makeVersionString( ${CULGT_VERSION} ${CULGT_VERSION_VARIABLE} )                
        
        # write new version file only if different from old version (avoids recompile if nothing changed)
        if( NOT OLD_VERSION_STRING STREQUAL NEW_VERSION_STRING )
            writeVersionString()
        endif()
    
    else()
        message( "No cuLGT version number available" )
        makeVersionString( "undefined" ${CULGT_VERSION_VARIABLE} )       
        writeVersionString() 
    endif()
    
else( IS_GIT_DIRECTORY )
    # if we are not git managed, a version file should come with the distribution (if not set version to undefined)
    if( NOT EXISTS ${PATH_TO_CURRENT_PROJECT}/${CULGT_VERSION_FILE} )
        message( "No cuLGT version number available" )
        makeVersionString( "undefined" ${CULGT_VERSION_VARIABLE} )       
        writeVersionString() 
    endif( NOT EXISTS ${PATH_TO_CURRENT_PROJECT}/${CULGT_VERSION_FILE} )
endif( IS_GIT_DIRECTORY )