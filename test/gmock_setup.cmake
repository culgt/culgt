#TODO probably should be findpackage style

#gmock
set( GMOCK_HOME $ENV{GMOCK_HOME} )
if( NOT DEFINED GMOCK_HOME )
    message( FATAL_ERROR "Please set GMOCK_HOME (in cmake or as environmental variable)" )
endif()

set( GMOCK_LIBDIR $ENV{GMOCK_LIBDIR} )
if( NOT DEFINED GMOCK_LIBDIR )
    if( EXISTS "${GMOCK_HOME}/build/libgmock.a")
        set( GMOCK_LIBDIR ${GMOCK_HOME}/build )
    else()
        message( FATAL_ERROR "Please set GMOCK_LIBDIR (in cmake or as environmental variable)" )
    endif()
endif()


#gtest
set( GTEST_HOME $ENV{GTEST_HOME} )
if( NOT DEFINED GTEST_HOME )
    message( FATAL_ERROR "Please set GTEST_HOME (in cmake or as environmental variable)" )
endif()

set( GTEST_LIBDIR $ENV{GTEST_LIBDIR} )
if( NOT DEFINED GTEST_LIBDIR )
    if( EXISTS "${GTEST_HOME}/build/libgtest.a")
        set( GTEST_LIBDIR ${GTEST_HOME}/build )
    else()
        message( FATAL_ERROR "Please set GTEST_LIBDIR (in cmake or as environmental variable)" )
    endif()
endif()