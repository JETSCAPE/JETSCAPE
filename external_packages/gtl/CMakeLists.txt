cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(GTL)

cmake_policy(SET CMP0042 NEW)

message("Compiler and Linker flags ...")

set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -O -fPIC -pipe -Wall -g -Wno-unused-variable -Wno-unused-private-field -Wno-reorder -Wno-logical-op-parentheses"
)

if(LINUX)
  message(STATUS "Linux : " ${CMAKE_HOST_SYSTEM})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

if(APPLE)
  message(STATUS "Apple : " ${CMAKE_HOST_SYSTEM})
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -Wno-builtin-macro-redefined -std=c++11 -stdlib=libc++"
    )
  elseif()
    set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -Wno-builtin-macro-redefined -std=c++11")
  endif()
  message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})
  message(STATUS "LD_FLAGS = " ${CMAKE_EXE_LINKER_FLAGS})
endif()

# w/o not compiling library automatically
add_subdirectory(./src)
link_directories(./lib)

# -----------------------------------------------------------
