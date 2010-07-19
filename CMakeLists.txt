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


# --- BOOST ---------------------------------------------------
SET(Boost_THREAD_FOUND 0)
FIND_PACKAGE(Boost COMPONENTS thread)
IF(Boost_THREAD_FOUND)

  # Look for the type of pattern to find.
  SET(BOOST_PATTERN_TOFIND "")
  IF(${CMAKE_BUILD_TYPE} MATCHES RELEASE)
    SET(BOOST_PATTERN_TOFIND optimized)
  ENDIF(${CMAKE_BUILD_TYPE} MATCHES RELEASE)
  IF(${CMAKE_BUILD_TYPE} MATCHES DEBUG)
    SET(BOOST_PATTERN_TOFIND debug)
  ENDIF(${CMAKE_BUILD_TYPE} MATCHES DEBUG)

  IF(BOOST_PATTERN_TOFIND)
    SET(BOOST_THREAD_LIST_LENGTH 0)
    LIST(LENGTH Boost_THREAD_LIBRARY BOOST_THREAD_LIST_LENGTH)
    # Debug or optimized library.
    IF(${BOOST_THREAD_LIST_LENGTH} MATCHES 4)
      LIST(FIND Boost_THREAD_LIBRARY ${BOOST_PATTERN_TOFIND} BTLB_index)
      LIST(REMOVE_AT Boost_THREAD_LIBRARY ${BTLB_index})
      LIST(GET Boost_THREAD_LIBRARY ${BTLB_index} BOOST_THREAD_LIB_NAME)
    # Default strategy
    ELSE(${BOOST_THREAD_LIST_LENGTH} MATCHES 4)
      SET(BOOST_THREAD_LIB_NAME ${Boost_THREAD_LIBRARY})
    ENDIF(${BOOST_THREAD_LIST_LENGTH} MATCHES 4)
  # Default strategy
  ELSE(BOOST_PATTERN_TOFIND)
    SET(BOOST_THREAD_LIB_NAME  ${Boost_THREAD_LIBRARY})
  ENDIF(BOOST_PATTERN_TOFIND)

ELSE(Boost_THREAD_FOUND)
  MESSAGE("WARNING: Boost thread not found. Trying without the thread library")
  FIND_PACKAGE(Boost)
  IF(NOT Boost_FOUND)
	MESSAGE(SEND_ERROR "!! ERROR: Boost not found.")
  ENDIF(NOT Boost_FOUND)
ENDIF(Boost_THREAD_FOUND)
MESSAGE(STATUS "Boost thread lib name: ${BOOST_THREAD_LIB_NAME} ${Boost_THREAD_LIBRARY}")

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

# --- LAPACK --------------------------------------------------
enable_language(Fortran)
SET(LAPACK_FOUND 0)
FIND_PACKAGE(LAPACK)
IF(LAPACK_FOUND)
#  IF (LINUX)
#    LIST(APPEND ${PROJECT_NAME}_LINK_FLAGS -lgfortran)
#  ENDIF(LINUX)
  LIST(APPEND ${PROJECT_NAME}_src_LDFLAGS ${LAPACK_LIBRARIES})
#  MESSAGE("LAPACK = "${${PROJECT_NAME}_src_LDFLAGS})
ENDIF(LAPACK_FOUND)

# --- BOOST NUMERIC BINDINGS-------------------------------------
IF(BOOSTNUMERICBINDINGS_DIR)
  SET(BOOSTNUMERICBINDINGS_DIR_SEARCH ${BOOSTNUMERICBINDINGS_DIR}/include)
ENDIF(BOOSTNUMERICBINDINGS_DIR)
FIND_PATH(BoostNumericBindings_INCLUDE_DIR
  NAMES boost/numeric/bindings/traits/config.hpp
  PATH_SUFFIXES ${SUFFIX_FOR_PATH}
  PATHS ${BOOSTNUMERICBINDINGS_DIR_SEARCH}
  DOC "The directory containing the boost bindings for lapack"
)
MESSAGE(STATUS "BoostNumericBindings_INCLUDE_DIR=${BoostNumericBindings_INCLUDE_DIR}")
IF(BoostNumericBindings_INCLUDE_DIR)
  SET(BoostNumericBindings_FOUND 1)
  INCLUDE_DIRECTORIES(${BoostNumericBindings_INCLUDE_DIR})
ENDIF(BoostNumericBindings_INCLUDE_DIR)


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
ADD_SUBDIRECTORY(unitTesting)

#
# Documentation produced by doxygen
#
INCLUDE(FindDoxygen)
IF(DOXYGEN_FOUND)
  add_subdirectory(doc)
ENDIF(DOXYGEN_FOUND)
