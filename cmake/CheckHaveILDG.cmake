# checks if required libraries are available

find_package(tinyxml QUIET )
include_directories( ${tinyxml_INCLUDE_DIRS} )

find_package( lime QUIET )
include_directories( ${lime_INCLUDE_DIRS} )

set( HAVE_ILDG_CHECKER TRUE )
set( ILDG_WARN_MESSAGE "ILDG file format disabled, missing packages: " )

if( NOT lime_FOUND )
    set( HAVE_ILDG_CHECKER FALSE )
    set( ILDG_WARN_MESSAGE ${ILDG_WARN_MESSAGE} "lime " )
endif()

if( NOT tinyxml_FOUND )
    set( HAVE_ILDG_CHECKER FALSE )
    set( ILDG_WARN_MESSAGE ${ILDG_WARN_MESSAGE} "tinyxml " )
endif()

if( HAVE_ILDG_CHECKER )
    message( STATUS "ILDG file format enabled" ) 
    set( CULGT_HAVE_LINKFILE_ILDG ON )
else()
    message( WARNING ${ILDG_WARN_MESSAGE} )
    set( CULGT_HAVE_LINKFILE_ILDG OFF )
endif()

configure_file( ${CULGT_HOME}/include/lattice/filetypes/filetype_config.h.in ${CULGT_HOME}/include/lattice/filetypes/filetype_config.h)