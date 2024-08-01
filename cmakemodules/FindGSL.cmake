#
# Try to find gnu scientific library GSL (see http://www.gnu.org/software/gsl/)
# Once run this will define:
#
# GSL_FOUND       = system has GSL lib
#
# GSL_LIBRARIES   = full path to the libraries on Unix/Linux with additional
# linker flags from "gsl-config --libs"
#
# CMAKE_GSL_CXX_FLAGS  = Unix compiler flags for GSL, essentially "`gsl-config
# --cxxflags`"
#
# GSL_INCLUDE_DIR      = where to find headers
#
# GSL_LINK_DIRECTORIES = link directories, useful for rpath on Unix
# GSL_EXE_LINKER_FLAGS = rpath on Unix
#
# Felix Woelk 07/2004 minor corrections Jan Woetzel
#
# www.mip.informatik.uni-kiel.de
# --------------------------------
#

if(WIN32)

  set(GSL_MINGW_PREFIX "c:/msys/local")
  set(GSL_MSVC_PREFIX "$ENV{LIB_DIR}")
  find_library(GSL_LIB gsl PATHS ${GSL_MINGW_PREFIX}/lib ${GSL_MSVC_PREFIX}/lib)
  # MSVC version of the lib is just called 'cblas'
  find_library(GSLCBLAS_LIB gslcblas cblas PATHS ${GSL_MINGW_PREFIX}/lib
                                                 ${GSL_MSVC_PREFIX}/lib)

  find_path(GSL_INCLUDE_DIR gsl/gsl_blas.h ${GSL_MINGW_PREFIX}/include
            ${GSL_MSVC_PREFIX}/include)

  if(GSL_LIB AND GSLCBLAS_LIB)
    set(GSL_LIBRARIES ${GSL_LIB} ${GSLCBLAS_LIB})
  endif(GSL_LIB AND GSLCBLAS_LIB)

else(WIN32)
  if(UNIX)
    set(GSL_CONFIG_PREFER_PATH
        "$ENV{GSL_HOME}/bin"
        CACHE STRING "preferred path to GSL (gsl-config)")
    find_program(GSL_CONFIG gsl-config ${GSL_CONFIG_PREFER_PATH}
                 $ENV{LIB_DIR}/bin /usr/local/bin/ /usr/bin/)
    # MESSAGE("DBG GSL_CONFIG ${GSL_CONFIG}")

    if(GSL_CONFIG)
      # set CXXFLAGS to be fed into CXX_FLAGS by the user:
      exec_program(
        ${GSL_CONFIG} ARGS
        --cflags
        OUTPUT_VARIABLE GSL_CXX_FLAGS)

      # set INCLUDE_DIRS to prefix+include
      exec_program(
        ${GSL_CONFIG} ARGS
        --prefix
        OUTPUT_VARIABLE GSL_PREFIX)
      set(GSL_INCLUDE_DIR
          ${GSL_PREFIX}/include
          CACHE STRING INTERNAL)

      # set link libraries and link flags
      exec_program(
        ${GSL_CONFIG} ARGS
        --libs
        OUTPUT_VARIABLE GSL_LIBRARIES)

      # extract link dirs for rpath
      exec_program(
        ${GSL_CONFIG} ARGS
        --libs
        OUTPUT_VARIABLE GSL_CONFIG_LIBS)

      # split off the link dirs (for rpath) use regular expression to match
      # wildcard equivalent "-L*<endchar>" with <endchar> is a space or a
      # semicolon
      string(REGEX MATCHALL "[-][L]([^ ;])+" GSL_LINK_DIRECTORIES_WITH_PREFIX
                   "${GSL_CONFIG_LIBS}")
      # MESSAGE("DBG
      # GSL_LINK_DIRECTORIES_WITH_PREFIX=${GSL_LINK_DIRECTORIES_WITH_PREFIX}")

      # remove prefix -L because we need the pure directory for LINK_DIRECTORIES

      if(GSL_LINK_DIRECTORIES_WITH_PREFIX)
        string(REGEX REPLACE "[-][L]" "" GSL_LINK_DIRECTORIES
                             ${GSL_LINK_DIRECTORIES_WITH_PREFIX})
      endif(GSL_LINK_DIRECTORIES_WITH_PREFIX)
      set(GSL_EXE_LINKER_FLAGS
          "-Wl,-rpath,${GSL_LINK_DIRECTORIES}"
          CACHE STRING INTERNAL)
      # MESSAGE("DBG  GSL_LINK_DIRECTORIES=${GSL_LINK_DIRECTORIES}")
      # MESSAGE("DBG  GSL_EXE_LINKER_FLAGS=${GSL_EXE_LINKER_FLAGS}")

      # ADD_DEFINITIONS("-DHAVE_GSL") SET(GSL_DEFINITIONS "-DHAVE_GSL")
      mark_as_advanced(GSL_CXX_FLAGS GSL_INCLUDE_DIR GSL_LIBRARIES
                       GSL_LINK_DIRECTORIES GSL_DEFINITIONS)

    else(GSL_CONFIG)

      if(GSL_FIND_REQUIRED)
        message(
          FATAL_ERROR
            "Could not find gsl-config. Please set it manually. GSL_CONFIG=${GSL_CONFIG}"
        )
      else(GSL_FIND_REQUIRED)
        message(STATUS "Could not find GSL")
        # TODO: Avoid cmake complaints if GSL is not found
      endif(GSL_FIND_REQUIRED)

    endif(GSL_CONFIG)

  endif(UNIX)
endif(WIN32)

if(GSL_LIBRARIES)
  if(GSL_INCLUDE_DIR OR GSL_CXX_FLAGS)

    set(GSL_FOUND 1)

    message(STATUS "Using GSL from ${GSL_PREFIX}")
    message(STATUS "  GSL_INCLUDE_DIR : ${GSL_INCLUDE_DIR}")
    message(STATUS "  GSL_CXX_FLAGS : ${GSL_CXX_FLAGS}")

  endif(GSL_INCLUDE_DIR OR GSL_CXX_FLAGS)
endif(GSL_LIBRARIES)
