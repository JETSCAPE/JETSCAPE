cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project (libJetScapeReader)

#for MacOSX
cmake_policy(SET CMP0042 NEW)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "../lib/")

set (SOURCES JetScapeReader.cc StringTokenizer.cc)

add_library(JetScapeReader SHARED ${SOURCES})

target_link_libraries(JetScapeReader JetScape JetScapeHydro GTL)


