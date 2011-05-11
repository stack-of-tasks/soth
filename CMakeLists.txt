#
#  Copyright 2010 CNRS
#

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/cpack.cmake)

SET(PROJECT_NAME soth)
SET(PROJECT_DESCRIPTION "...")
SET(PROJECT_URL "")

# Disable -Werror on Unix for now.
SET(CXX_DISABLE_WERROR True)

project(${PROJECT_NAME})

SET(PROJECT_VERSION 0.0.1)

SET(CMAKE_VERBOSE_MAKEFILE true)

SETUP_PROJECT()


#
# Look for dependencies
#
SET(CMAKE_MODULE_PATH ${${PROJECT_NAME}_SOURCE_DIR}/CMakeModules)

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
set(EIGEN3_INCLUDE_DIR "C:/Users/Keith/Libraries/Eigen/include/eigen3")
FIND_PACKAGE(Eigen3)
IF(EIGEN3_FOUND)
  INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
  PKG_CONFIG_APPEND_CFLAGS(${EIGEN3_INCLUDE_DIR} -DEIGEN2_SUPPORT)
  IF(WIN32)
	PKG_CONFIG_APPEND_CFLAGS(${EIGEN3_INCLUDE_DIR} "-DNOMINMAX")
  ENDIF(WIN32)

ELSE()
  MESSAGE(
    FATAL_ERROR "Failed to detect Eigen3\n"
  )
ENDIF(EIGEN3_FOUND)

# --- BOOST HEADER ---------------------------------------------------
# Search for dependencies.
# Boost
SET(BOOST_COMPONENTS signals)
SEARCH_FOR_BOOST()

# ADD_SUBDIRECTORY(doc)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(unitTesting)

SETUP_PROJECT_FINALIZE()
SETUP_PROJECT_CPACK()
