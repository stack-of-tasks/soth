#
#  Copyright 2010 CNRS
#
#  Author: Florent Lamiraux
#

SET(LIBRARY_NAME ${PROJECT_NAME})


SET(${PROJECT_NAME}_SRC
    soth.cpp
    Algebra.cpp
    Bound.cpp
    Stage.cpp
    BaseY.cpp
    HCOD.cpp
    Givens.cpp
    debug.cpp
)


# ----------------------------------------------------
# --- LIBS -------------------------------------------
# ----------------------------------------------------

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)

ADD_LIBRARY(${LIBRARY_NAME}
  SHARED
  ${${PROJECT_NAME}_SRC})

SET_TARGET_PROPERTIES(${LIBRARY_NAME}
  PROPERTIES
  SOVERSION ${PROJECT_VERSION}
  INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)

LIST2STRING(_cf "${${PROJECT_NAME}_src_CXXFLAGS}")
LIST2STRING(_lf "${${PROJECT_NAME}_src_LDFLAGS}")
SET_TARGET_PROPERTIES(${PROJECT_NAME_LIB}
  PROPERTIES
  COMPILE_FLAGS "${_cf}"
  LINK_FLAGS "${_lf}"
)

INSTALL(TARGETS ${LIBRARY_NAME}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
