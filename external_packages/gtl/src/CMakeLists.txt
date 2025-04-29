cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project (GTL_LIB)

#for MacOSX
cmake_policy(SET CMP0042 NEW)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "../lib")

include_directories(../include)

#set (SOURCES "*.cpp")

file(GLOB SOURCES "./*.cpp")

add_library(GTL SHARED ${SOURCES})
