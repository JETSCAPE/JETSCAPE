# Install script for directory: /home/kevin/DukeQCD/JETSCAPE/external_packages/LBTD/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/kevin/DukeQCD/JETSCAPE/build")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/libLBTD" TYPE STATIC_LIBRARY FILES "/home/kevin/DukeQCD/JETSCAPE/build/external_packages/LBTD/src/libLBTD.a")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/libLBTD" TYPE FILE FILES
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/LBTD/src/logo.h"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/LBTD/src/lorentz.h"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/LBTD/src/stat.h"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/LBTD/src/Rate.h"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/LBTD/src/matrix_elements.h"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/LBTD/src/workflow.h"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/LBTD/src/vwrapper.h"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/LBTD/src/sampler.h"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/LBTD/src/TableBase.h"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/LBTD/src/converged.h"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/LBTD/src/minimizer.h"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/LBTD/src/cubature.h"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/LBTD/src/StochasticBase.h"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/LBTD/src/predefine.h"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/LBTD/src/simpleLogger.h"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/LBTD/src/random.h"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/LBTD/src/approx_functions.h"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/LBTD/src/fast_exp.h"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/LBTD/src/integrator.h"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/LBTD/src/Srandom.h"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/LBTD/src/Xsection.h"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/LBTD/src/Langevin.h"
    )
endif()

