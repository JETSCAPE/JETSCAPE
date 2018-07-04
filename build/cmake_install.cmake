# Install script for directory: /home/kevin/DukeQCD/JETSCAPE

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
    set(CMAKE_INSTALL_CONFIG_NAME "")
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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY FILES
    "/home/kevin/DukeQCD/JETSCAPE/src/framework/"
    "/home/kevin/DukeQCD/JETSCAPE/src/hadronization/"
    "/home/kevin/DukeQCD/JETSCAPE/src/initialstate/"
    "/home/kevin/DukeQCD/JETSCAPE/src/hydro/"
    "/home/kevin/DukeQCD/JETSCAPE/src/jet/"
    "/home/kevin/DukeQCD/JETSCAPE/src/reader/"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/iSS/src/"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/hydro_from_external_file/src/"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/music/src/"
    "/home/kevin/DukeQCD/JETSCAPE/external_packages/trento/src/"
    FILES_MATCHING REGEX "/[^/]*\\.h[^/]*$" REGEX "/gtl$" EXCLUDE REGEX "/trento$" EXCLUDE REGEX "/tests$" EXCLUDE REGEX "/data\\_table$" EXCLUDE REGEX "/LBT\\-tables$" EXCLUDE REGEX "/Martini$" EXCLUDE REGEX "/googletest$" EXCLUDE)
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY FILES "/home/kevin/DukeQCD/JETSCAPE/external_packages/gtl/include/" FILES_MATCHING REGEX "/GTL\\/[^/]*\\.h[^/]*$")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE DIRECTORY FILES
    "/home/kevin/DukeQCD/JETSCAPE/build/src/lib/"
    "/home/kevin/DukeQCD/JETSCAPE/build/lib/"
    "/home/kevin/DukeQCD/JETSCAPE/build/external_packages/gtl/lib/"
    "/home/kevin/DukeQCD/JETSCAPE/build/external_packages/iSS/src/"
    "/home/kevin/DukeQCD/JETSCAPE/build/external_packages/music/src/"
    "/home/kevin/DukeQCD/JETSCAPE/build/external_packages/trento/src/"
    FILES_MATCHING REGEX "/lib[^/]*\\.[^/]*$" REGEX "/CMakeFiles$" EXCLUDE)
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/kevin/DukeQCD/JETSCAPE/build/external_packages/trento/cmake_install.cmake")
  include("/home/kevin/DukeQCD/JETSCAPE/build/external_packages/LBTD/cmake_install.cmake")
  include("/home/kevin/DukeQCD/JETSCAPE/build/external_packages/cmake_install.cmake")
  include("/home/kevin/DukeQCD/JETSCAPE/build/external_packages/gtl/cmake_install.cmake")
  include("/home/kevin/DukeQCD/JETSCAPE/build/src/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/kevin/DukeQCD/JETSCAPE/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
