cmake_minimum_required(VERSION 2.8)
project(lib_clvisc)

#for MacOSX
cmake_policy(SET CMP0042 NEW)

if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external_packages/PyVisc")
    message(FATAL_ERROR "Error: CLVisc source has not been downloaded in external_packages by ./get_clvisc.sh")
else()
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/eos_table/s95_pce165.dat DESTINATION ${CMAKE_BINARY_DIR}/data_table)
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../PyVisc/pyvisc/kernel/ DESTINATION ${CMAKE_BINARY_DIR}/clvisc_kernel)
endif()

set ( SRC_LIST 
    ${PROJECT_SOURCE_DIR}/include/cl.hpp
    ${PROJECT_SOURCE_DIR}/opencl_backend.cc
    ${PROJECT_SOURCE_DIR}/include/clideal.h
    ${PROJECT_SOURCE_DIR}/clideal.cc
    ${PROJECT_SOURCE_DIR}/include/clvisc.h
    ${PROJECT_SOURCE_DIR}/clvisc.cc
    ${PROJECT_SOURCE_DIR}/test_visc.cc
    #${PROJECT_SOURCE_DIR}/test.cc
    )
set( CMAKE_CXX_COMPILER  "g++" )

# will enable more verbose opencl error message 
# instead of one int error id
set(ENV{CL_LOG_ERRORS} "stderr")

set ( CMAKE_USE_RELATIVE_PATH "true")
set ( CMAKE_EXECUTABLE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build/ )
set ( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build/ )

set ( BSZ 256 )   

set ( DUMPDATA_FROM_GPU "true" )

option ( USE_DEVICE_GPU "Use GPU parallel" on )
option ( USE_SINGLE_PRECISION "Use float instead of double" on )
option ( CL_USE_DEPRECATED_OPENCL_1_1_APIS "Use Nvidia GPU" off )

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
        message( "No OpenCL CPP bindings found, will use clvisc_wrapper/include/cl.hpp" )
    endif( OPENCL_HAS_CPP_BINDINGS )
    include_directories( ${OPENCL_INCLUDE_DIRS} )
    set ( EXTRA_LIBS ${EXTRA_LIBS} ${OPENCL_LIBRARIES} )
else( OPENCL_FOUND )
    #### specific for GSI green cube
    message( "CMake did not find OpenCL, set by hand" )
endif( OPENCL_FOUND )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wno-deprecated")
include_directories( ${PROJECT_SOURCE_DIR}/include )

### add_executable ( test_clvisc ${SRC_LIST} )
### target_link_libraries ( test_clvisc ${EXTRA_LIBS} )

#### compile the clvisc to one library that will be used by clvisc+jetscape interface
add_library(lib_clvisc SHARED ${SRC_LIST})
target_link_libraries (lib_clvisc ${EXTRA_LIBS})
