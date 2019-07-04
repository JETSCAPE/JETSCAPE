cmake_minimum_required(VERSION 2.8)
project(Lido CXX)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" CXX11)
if(NOT CXX11)
  message(FATAL_ERROR "Your compiler does not support C++11")
endif()

# Find and use Boost.
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
find_package(Boost 1.5.0 REQUIRED COMPONENTS filesystem program_options system log log_setup thread REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

# Search for the HDF5 C++ library.  The CMake FindHDF5 package does not work
# correctly for C++, because it will report success if it finds only the C
# library and not C++.  It also uses the rather unfortunate 'h5c++' compiler
# wrapper which sets a lot of unnecessary crap.

find_library(HDF5_CXX_LIB hdf5_cpp)
if("${HDF5_CXX_LIB}" STREQUAL "HDF5_CXX_LIB-NOTFOUND")
  message(STATUS "HDF5 C++ library not found")
else()
  message(STATUS "Found HDF5 C++ library: ${HDF5_CXX_LIB}")
  # Having found the C++ lib, also locate the header and C lib.
  find_path(HDF5_INCLUDE_DIRS H5Cpp.h)
  include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})
  find_library(HDF5_C_LIB hdf5)
  message(STATUS "Found HDF5 C library: ${HDF5_C_LIB}")
  # Also link to zlib for HDF5 compression.
  set(HDF5_LIBRARIES ${HDF5_CXX_LIB} ${HDF5_C_LIB} z)
  # Set preprocessor definition so HDF5 output is compiled.
  add_definitions(-DTRENTO_HDF5)
endif()

# find and use gsl
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake)
FIND_PACKAGE(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})


set(PROJECT_VERSION_MAJOR 1)
set(PROJECT_VERSION_MINOR 1)
set(PROJECT_VERSION_PATCH 1)
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

# Set default architecture flags to optimize for the current system,
# but don't override user settings.
if(NOT CMAKE_CXX_FLAGS MATCHES "-march|-mtune")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

# warnings adapted from http://stackoverflow.com/a/9862800
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Werror")

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og")
endif()

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  # Check for sanitizer support.
  # (see e.g. https://code.google.com/p/address-sanitizer/wiki/AddressSanitizer)
  # check_cxx_compiler_flag() doesn't work here because it doesn't set the test
  # flag at link time, so use check_cxx_source_compiles() instead.
  include(CheckCXXSourceCompiles)
  set(CMAKE_REQUIRED_QUIET TRUE)

  # First try both AddressSanitizer and UndefinedBehaviorSanitizer.
  set(CMAKE_REQUIRED_FLAGS "-fsanitize=address,undefined")
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

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_REQUIRED_FLAGS} --coverage")

set(LIBRARY_NAME "lib${PROJECT_NAME}")
add_subdirectory(src)
# add_subdirectory(test)
# add_subdirectory(doc)
