#
#  Copyright 2010 CNRS
#

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/cpack.cmake)
INCLUDE(cmake/eigen.cmake)

SET(PROJECT_NAME soth)
SET(PROJECT_DESCRIPTION "...")
SET(PROJECT_URL "http://github.com/laas/soth")

# Disable -Werror on Unix for now.
SET(CXX_DISABLE_WERROR True)
SET(CMAKE_VERBOSE_MAKEFILE true)

project(${PROJECT_NAME})

SETUP_PROJECT()
SEARCH_FOR_EIGEN()

# Handle OS specificities 
INCLUDE(CheckIncludeFiles)
CHECK_INCLUDE_FILES("sys/time.h" SYS_TIME_H)
IF(SYS_TIME_H)
ADD_DEFINITIONS("-DHAVE_SYS_TIME_H")
ENDIF(SYS_TIME_H)

# --- BOOST HEADER ---------------------------------------------------
# Search for dependencies.
# Boost
SET(BOOST_COMPONENTS program_options signals)
SEARCH_FOR_BOOST()

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(unitTesting)

PKG_CONFIG_APPEND_CFLAGS(${_Eigen_CFLAGS})
PKG_CONFIG_APPEND_LIBS("soth")

SETUP_PROJECT_FINALIZE()
SETUP_PROJECT_CPACK()

