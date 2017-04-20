##############
# Unit Tests
##############

SET(jetclass_src JetClass.cpp)
# compile JetClass to library to be used in other modules
add_library(jetclass_lib STATIC ${jetclass_src})

if (test)
  add_subdirectory(test)
endif()
