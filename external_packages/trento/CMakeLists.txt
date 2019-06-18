
cmake_minimum_required(VERSION 3.4)
project(trento VERSION 1.5.1 LANGUAGES CXX)



# require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(STATIC "build static executable" OFF)
if(CMAKE_CXX_FLAGS MATCHES "(^| )-static($| )")
  set(STATIC ON)
endif()
# adapted from cmake source: Modules/Compiler/CrayPrgEnv.cmake
if(STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(CMAKE_LINK_SEARCH_START_STATIC TRUE)
endif()


#old
#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#include(CheckCXXCompilerFlag)
#check_cxx_compiler_flag("-std=c++11" CXX11)
#if(NOT CXX11)
#message(FATAL_ERROR "Your compiler does not support C++11")
#endif()


# Find and use Boost.
find_package(Boost 1.60 REQUIRED COMPONENTS filesystem program_options system)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
message(STATUS "${Boost_INCLUDE_DIRS}, ${Boost_LIBRARIES}")
# Find and use GSL.
#set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/src/CMake) 
FIND_PACKAGE(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIR} ${GSLCBLAS_INCLUDE_DIR})
message(STATUS ${GSL_INCLUDE_DIR} ${GSLCBLAS_INCLUDE_DIR})
# Find the HDF5 C++ library.  The CMake FindHDF5 module does not work reliably,
# e.g. sometimes it links unnecessary libraries, other times it fails to include
# libz.  The following works on all systems I have tested (Ubuntu, Fedora, Arch,
# Cray Linux at NERSC, RHEL on OSG).

# Read HDF5_ROOT either as a CMake option or environment variable.
if(NOT HDF5_ROOT)
  set(HDF5_ROOT $ENV{HDF5_ROOT})
endif()

# If HDF5_ROOT was set, search for the header and libraries in that location,
# otherwise try to extract the locations from the h5c++ compiler wrapper.
if(HDF5_ROOT)
  set(HDF5_INCLUDE_HINT "${HDF5_ROOT}/include")
  set(HDF5_LIBRARY_HINT "${HDF5_ROOT}/lib")
else()
  find_program(H5CXX h5c++)
  if(H5CXX)
    if(STATIC)
      set(H5CXX_ARG "-noshlib")
    else()
      set(H5CXX_ARG "-shlib")
    endif()
    # parse output to extract include and library dirs (-I and -L flags)
    execute_process(COMMAND ${H5CXX} -show ${H5CXX_ARG} OUTPUT_VARIABLE H5CXX_OUTPUT)
    foreach(TYPE "INCLUDE" "LIBRARY")
      string(SUBSTRING ${TYPE} 0 1 FLAG)
      set(FLAG " -${FLAG}")
      string(REGEX MATCH "${FLAG}[^ ]+" FLAG_HINT ${H5CXX_OUTPUT})
      string(REPLACE "${FLAG}" "" HDF5_${TYPE}_HINT "${FLAG_HINT}")
    endforeach()
  endif()
endif()

# Find all required HDF5 components.
find_path(HDF5_INCLUDE_DIRS H5Cpp.h HINTS ${HDF5_INCLUDE_HINT})
find_library(HDF5_CXX_LIB hdf5_cpp HINTS ${HDF5_LIBRARY_HINT})
find_library(HDF5_C_LIB hdf5 HINTS ${HDF5_LIBRARY_HINT})
find_library(ZLIB z)

# Group HDF5 components by whether they were found or not.
set(HDF5_VARS_FOUND "")
set(HDF5_VARS_MISSING "")
set(HDF5_FOUND TRUE)
foreach(VAR HDF5_INCLUDE_DIRS HDF5_CXX_LIB HDF5_C_LIB ZLIB)
  if(${VAR})
    set(HDF5_VARS_FOUND "${HDF5_VARS_FOUND} ${${VAR}}")
  else()
    set(HDF5_FOUND FALSE)
    set(HDF5_VARS_MISSING "${HDF5_VARS_MISSING} ${VAR}")
  endif()
endforeach()

# Report HDF5 status and set relevant variables if found.
if(HDF5_FOUND)
  set(HDF5_LIBRARIES ${HDF5_CXX_LIB} ${HDF5_C_LIB} ${ZLIB})
  include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})
  add_definitions(-DTRENTO_HDF5)
  message(STATUS "Found HDF5:${HDF5_VARS_FOUND}")
else()
  message(STATUS "HDF5 not found (missing:${HDF5_VARS_MISSING}, found:${HDF5_VARS_FOUND})")
endif()


set(PROJECT_VERSION_MAJOR 1)
set(PROJECT_VERSION_MINOR 0)
set(PROJECT_VERSION_PATCH 0)
set(PROJECT_VERSION
  "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")

# default build type: Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# default install prefix: ~/.local
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local"
      CACHE PATH "Install path prefix, prepended onto install directories."
      FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")

# By default, optimize for the system's native architecture.  Disable via the
# NATIVE option.  In addition, detect if another architecture flag is already
# set and do not override it.
option(NATIVE "compile for native architecture" ON)
if(NATIVE)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    # Intel compiler: search for -m, -x, -ax flags; if not found add -xHost.
    if(NOT CMAKE_CXX_FLAGS MATCHES "(^| )-(m|x|ax)[^ ]+")
      string(APPEND CMAKE_CXX_FLAGS " -xHost")
    endif()
  else()
    # All other compilers (gcc-like): search for -m flags;
    # if not found add -march=native.
    if(NOT CMAKE_CXX_FLAGS MATCHES "(^| )-m[^ ]+")
      string(APPEND CMAKE_CXX_FLAGS " -march=native")
    endif()
  endif()
endif()

string(APPEND CMAKE_CXX_FLAGS " -Wall -Wextra")

option(MORE_WARNINGS "enable more compiler warnings" ON)
if(MORE_WARNINGS AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  # adapted from http://stackoverflow.com/a/9862800
  string(APPEND CMAKE_CXX_FLAGS " -Wpedantic -Wcast-align -Wcast-qual -Wdisabled-optimization -Wformat=2 -Wmissing-declarations -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wsign-conversion -Wsign-promo -Wstrict-overflow=2 -Wundef")
endif()

# disable silly clang warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  string(APPEND CMAKE_CXX_FLAGS " -Wno-missing-braces -Wno-c++11-narrowing")
endif()

string(APPEND CMAKE_CXX_FLAGS_DEBUG " -Werror")

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  string(APPEND CMAKE_CXX_FLAGS_DEBUG " -Og")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  # Check for sanitizer support.
  # (see e.g. https://code.google.com/p/address-sanitizer/wiki/AddressSanitizer)
  # check_cxx_compiler_flag() doesn't work here because it doesn't set the test
  # flag at link time, so use check_cxx_source_compiles() instead.
  include(CheckCXXSourceCompiles)
  set(CMAKE_REQUIRED_QUIET TRUE)

  # First try both AddressSanitizer and UndefinedBehaviorSanitizer.
  # Disable the vptr sanitizer because it throws errors with boost::any.
  set(CMAKE_REQUIRED_FLAGS "-fsanitize=address,undefined -fno-sanitize=vptr")
  check_cxx_source_compiles("int main() { return 0; }" AddrUndefSanitizer)

  # Lacking that, try AddressSanitizer only.
  if(NOT AddrUndefSanitizer)
    set(CMAKE_REQUIRED_FLAGS "-fsanitize=address")
    check_cxx_source_compiles("int main() { return 0; }" AddrSanitizer)
    if(NOT AddrSanitizer)
      # out of luck...
      set(CMAKE_REQUIRED_FLAGS "")
    endif()
  endif()
endif()

string(APPEND CMAKE_CXX_FLAGS_DEBUG " ${CMAKE_REQUIRED_FLAGS} --coverage")

set(LIBRARY_NAME "lib${PROJECT_NAME}")

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(doc)
