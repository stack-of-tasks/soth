#
#  Copyright 2010 CNRS
#  Author: Florent Lamiraux
#

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(PROJECT_NAME soth)
project(${PROJECT_NAME})

SET(PROJECT_VERSION 1.0)

SET(CMAKE_VERBOSE_MAKEFILE true)

IF(WIN32)
  SET(LIBDIR_KW "/LIBPATH:")
  SET(LIBINCL_KW "")
  SET(LIB_EXT ".lib")
ENDIF(WIN32)

IF(UNIX)
  SET(LIBDIR_KW "-L")
  SET(LIBINCL_KW "-l")
  SET(LIB_EXT "")
ENDIF(UNIX)

#
# Look for dependencies
#
SET(CMAKE_MODULE_PATH ${${PROJECT_NAME}_SOURCE_DIR}/CMakeModules)

# Load Pkgconfig macros
INCLUDE(FindPkgConfig)
SET(PACKAGE_REQUIREMENTS "")

# --- UTILS ---------------------------------------------------------
MACRO(LIST2STRING _VAR_ARG _LIST_ARG)
  SET(${_VAR_ARG})
  SET(_FIRST_STEP 1)
  FOREACH(_ITEM ${ARGV})
     IF(_FIRST_STEP)
	SET(_FIRST_STEP 0)
     ELSE(_FIRST_STEP)
           SET (${_VAR_ARG} "${${_VAR_ARG}} ${_ITEM}")
     ENDIF(_FIRST_STEP)
  ENDFOREACH(_ITEM)
ENDMACRO(LIST2STRING)


# --- EIGEN ----------------------------------------------------------
FIND_PACKAGE(Eigen3)
IF(EIGEN3_FOUND)
  INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
ENDIF(EIGEN3_FOUND)

# --- BOOST HEADER ---------------------------------------------------
FIND_PACKAGE( Boost 1.36.0 )
IF(Boost_FOUND)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
ENDIF()

#
# pkg-config file
#

SET(install_pkg_prefix "\${prefix}")
SET(install_pkg_exec_prefix "\${exec_prefix}")
SET(install_pkg_libdir "\${libdir}")
SET(install_pkg_datarootdir "\${datarootdir}")
SET(install_pkg_include_dir "\${includedir}")

IF(UNIX)
  SET (${PROJECT_NAME}_LDFLAGS "-Wl,-R${install_pkg_libdir} ${LIBINCL_KW}${PROJECT_NAME}${LIB_EXT} ${${PROJECT_NAME}_LDFLAGS}")
ENDIF(UNIX)

IF(WIN32)
  SET (${PROJECT_NAME}_LDFLAGS "${LIBINCL_KW}${PROJECT_NAME}${LIB_EXT} ${${PROJECT_NAME}_LDFLAGS}")
ENDIF(WIN32)

# Install the file in the appropriate place.
SET(install_libdir_pkgconfig ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}.pc.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc)

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc
  DESTINATION ${install_libdir_pkgconfig}
  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ OWNER_WRITE)

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(include)
ADD_SUBDIRECTORY(unitTesting)

#
# Documentation produced by doxygen
#
INCLUDE(FindDoxygen)
IF(DOXYGEN_FOUND)
  add_subdirectory(doc)
ENDIF(DOXYGEN_FOUND)
