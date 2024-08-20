# * Try to find OpenCL This module tries to find an OpenCL implementation on
#   your system. It supports AMD / ATI, Apple and NVIDIA implementations, but
#   should work, too.
#
# To set manually the paths, define these environment variables: OpenCL_INCPATH
# - Include path (e.g. OpenCL_INCPATH=/opt/cuda/4.0/cuda/include) OpenCL_LIBPATH
# - Library path (e.h. OpenCL_LIBPATH=/usr/lib64/nvidia)
#
# Once done this will define OPENCL_FOUND        - system has OpenCL
# OPENCL_INCLUDE_DIRS  - the OpenCL include directory OPENCL_LIBRARIES    - link
# these to use OpenCL
#
# WIN32 should work, but is untested

find_package(PackageHandleStandardArgs)

set(OPENCL_VERSION_STRING "0.1.0")
set(OPENCL_VERSION_MAJOR 0)
set(OPENCL_VERSION_MINOR 1)
set(OPENCL_VERSION_PATCH 0)

if(APPLE)

  find_library(OPENCL_LIBRARIES OpenCL DOC "OpenCL lib for OSX")
  find_path(OPENCL_INCLUDE_DIRS OpenCL/cl.h DOC "Include for OpenCL on OSX")
  find_path(_OPENCL_CPP_INCLUDE_DIRS OpenCL/cl.hpp
            DOC "Include for OpenCL CPP bindings on OSX")

else(APPLE)

  if(WIN32)

    find_path(OPENCL_INCLUDE_DIRS CL/cl.h)
    find_path(_OPENCL_CPP_INCLUDE_DIRS CL/cl.hpp)

    # The AMD SDK currently installs both x86 and x86_64 libraries This is only
    # a hack to find out architecture
    if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "AMD64")
      set(OPENCL_LIB_DIR "$ENV{ATISTREAMSDKROOT}/lib/x86_64")
    else(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "AMD64")
      set(OPENCL_LIB_DIR "$ENV{ATISTREAMSDKROOT}/lib/x86")
    endif(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "AMD64")
    find_library(OPENCL_LIBRARIES OpenCL.lib PATHS ${OPENCL_LIB_DIR} ENV
                                                   OpenCL_LIBPATH)

    get_filename_component(_OPENCL_INC_CAND ${OPENCL_LIB_DIR}/../../include
                           ABSOLUTE)

    # On Win32 search relative to the library
    find_path(OPENCL_INCLUDE_DIRS CL/cl.h PATHS "${_OPENCL_INC_CAND}" ENV
                                                OpenCL_INCPATH)
    find_path(_OPENCL_CPP_INCLUDE_DIRS CL/cl.hpp PATHS "${_OPENCL_INC_CAND}"
                                                       ENV OpenCL_INCPATH)

  else(WIN32)

    # Unix style platforms
    find_library(OPENCL_LIBRARIES OpenCL PATHS ENV LD_LIBRARY_PATH ENV
                                               OpenCL_LIBPATH)

    get_filename_component(OPENCL_LIB_DIR ${OPENCL_LIBRARIES} PATH)
    get_filename_component(_OPENCL_INC_CAND ${OPENCL_LIB_DIR}/../../include
                           ABSOLUTE)

    # The AMD SDK currently does not place its headers in /usr/include,
    # therefore also search relative to the library
    find_path(OPENCL_INCLUDE_DIRS CL/cl.h
              PATHS ${_OPENCL_INC_CAND} "/usr/local/cuda/include"
                    "/opt/AMDAPP/include" ENV OpenCL_INCPATH)
    find_path(_OPENCL_CPP_INCLUDE_DIRS CL/cl.hpp
              PATHS ${_OPENCL_INC_CAND} "/usr/local/cuda/include"
                    "/opt/AMDAPP/include" ENV OpenCL_INCPATH)

  endif(WIN32)

endif(APPLE)

find_package_handle_standard_args(OpenCL DEFAULT_MSG OPENCL_LIBRARIES
                                  OPENCL_INCLUDE_DIRS)

if(_OPENCL_CPP_INCLUDE_DIRS)
  set(OPENCL_HAS_CPP_BINDINGS TRUE)
  list(APPEND OPENCL_INCLUDE_DIRS ${_OPENCL_CPP_INCLUDE_DIRS})
  # This is often the same, so clean up
  list(REMOVE_DUPLICATES OPENCL_INCLUDE_DIRS)
endif(_OPENCL_CPP_INCLUDE_DIRS)

mark_as_advanced(OPENCL_INCLUDE_DIRS)
