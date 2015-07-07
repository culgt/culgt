# this module defines
# lime_FOUND
# lime_LIBRARIES
# lime_INCLUDE_DIRS

include( LibFindMacros )

find_path(lime_INCLUDE_DIR
  NAMES lime.h
)

find_library(lime_LIBRARY
  NAMES lime
)

libfind_process( lime )