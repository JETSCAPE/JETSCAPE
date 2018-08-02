project(clvisc_for_jetscape)

cmake_minimum_required(VERSION 2.8)

#cl_spec.cpp  EosPCEv0.cl  helper.h  kernel_ideal.cl
set ( SRC_LIST 
    ${PROJECT_SOURCE_DIR}/include/cl.hpp
    ${PROJECT_SOURCE_DIR}/opencl_backend.cc
    ${PROJECT_SOURCE_DIR}/include/clideal.h
    ${PROJECT_SOURCE_DIR}/clideal.cc
    ${PROJECT_SOURCE_DIR}/test.cc
    )
#aux_source_directory(src/ SRC_LIST)
#aux_source_directory(inc/ SRC_LIST)
set( CMAKE_CXX_COMPILER  "g++" )

set ( CMAKE_USE_RELATIVE_PATH "true")
set ( CMAKE_EXECUTABLE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build/ )
set ( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build/ )

set ( BSZ 256 )   

set ( DUMPDATA_FROM_GPU "true" )

option ( USE_DEVICE_GPU "Use GPU parallel" on )
option ( USE_SINGLE_PRECISION "Use float instead of double" on )
option ( CL_USE_DEPRECATED_OPENCL_1_1_APIS "Use Nvidia GPU" off )
option ( LOEWE_CSC "Use LOEWE AMD GPUs for spec kernel" on )


configure_file (
   ${PROJECT_SOURCE_DIR}/include/Config.h.in
   ${PROJECT_SOURCE_DIR}/include/Config.h
)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../../cmakemodules/")
 
find_package( OpenCL REQUIRED )
 
include_directories( "/usr/local/include/" )
link_directories( "/usr/local/lib" )

if( OPENCL_FOUND )
    if( OPENCL_HAS_CPP_BINDINGS )
        message( "OpenCL has CPP bindings. Full include is: " ${OPENCL_INCLUDE_DIRS} )
    else( OPENCL_HAS_CPP_BINDINGS )
        message( "No OpenCL CPP bindings found, will use inc/cl.hpp" )
    endif( OPENCL_HAS_CPP_BINDINGS )
    include_directories( ${OPENCL_INCLUDE_DIRS} )
    set ( EXTRA_LIBS ${EXTRA_LIBS} ${OPENCL_LIBRARIES} )
else( OPENCL_FOUND )
    #### specific for GSI green cube
    message( "CMake did not find OpenCL, set by hand" )
endif( OPENCL_FOUND )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

include_directories( ${PROJECT_SOURCE_DIR}/include )

#include_directories( ${PROJECT_SOURCE_DIR}/../ReadConfig/inc )
#LINK_DIRECTORIES( ${PROJECT_SOURCE_DIR}/../ReadConfig/lib )

add_executable ( test ${SRC_LIST} )
target_link_libraries ( test ${EXTRA_LIBS} )
