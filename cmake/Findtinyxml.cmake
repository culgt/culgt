# this module defines
# tinyxml_FOUND
# tinyxml_LIBRARIES
# tinyxml_INCLUDE_DIRS

include( LibFindMacros )

find_path(tinyxml_INCLUDE_DIR
  NAMES tinyxml.h
)

find_library(tinyxml_LIBRARY
  NAMES tinyxml
)

libfind_process( tinyxml )