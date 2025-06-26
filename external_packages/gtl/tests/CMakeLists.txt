cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project (glt_test)

cmake_policy(SET CMP0042 NEW)
cmake_policy(SET CMP0015 NEW)

message("Include Directories ...")

include_directories(../include)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "./")

message("Compiler and Linker flags ...")

set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -O -fPIC -pipe -Wall -g -Wno-unused-variable -Wno-unused-private-field -Wno-reorder -Wno-comment"
)

if(LINUX)
  message(STATUS "Linux : " ${CMAKE_HOST_SYSTEM})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

if(APPLE)
  message(STATUS "Apple : " ${CMAKE_HOST_SYSTEM})
  set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -Wno-macro-redefined -std=c++11 -stdlib=libc++")
  message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})
  message(STATUS "LD_FLAGS = " ${CMAKE_EXE_LINKER_FLAGS})
endif()

link_directories(../lib)

add_executable(my_test my_test.cc)
add_executable(my_test2 my_test2.cc)
add_executable(my_test3 my_test3.cc)
add_executable(my_test4 my_test4.cc)

target_link_libraries(my_test libGTL.dylib)
target_link_libraries(my_test2 libGTL.dylib)
target_link_libraries(my_test3 libGTL.dylib)
target_link_libraries(my_test4 libGTL.dylib)
