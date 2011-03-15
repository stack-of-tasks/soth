#
#  Copyright 2010 CNRS
#
#  Author: Florent Lamiraux
#

SET(LIBRARY_NAME ${PROJECT_NAME})


SET(${PROJECT_NAME}_SRC
    ActiveSet.cpp
    Allocator.cpp
    Algebra.cpp
    Bound.cpp
    BaseY.cpp
    HCOD.cpp
    Givens.cpp
    debug.cpp
    Stage.cpp
    BasicStage.cpp
)

SET(${PROJECT_NAME}_HEADERS
     ActiveSet.hpp
     Algebra.hpp
     Allocator.hpp
     BaseY.hpp
     BasicStage.hpp
     Bound.hpp
     debug.hpp
     DestructiveColPivQR.hpp
     Givens.hpp
     HCOD.hpp
     COD.hpp
     solvers.hpp
     Stage.hpp
     SubMatrix.hpp
)

# ----------------------------------------------------
# --- INCLUDE ----------------------------------------
# ----------------------------------------------------

MAKE_DIRECTORY("${${PROJECT_NAME}_BINARY_DIR}/include/soth")

FOREACH(header ${${PROJECT_NAME}_HEADERS})
  GET_FILENAME_COMPONENT(headerName ${header} NAME)
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                  ${${PROJECT_NAME}_SOURCE_DIR}/src/${header}
		  ${${PROJECT_NAME}_BINARY_DIR}/include/${PROJECT_NAME}/${headerName})
  INSTALL(FILES ${${PROJECT_NAME}_SOURCE_DIR}/${header}
	DESTINATION ${CMAKE_INSTALL_PREFIX}/include/${PROJECT_NAME}
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ OWNER_WRITE
)

ENDFOREACH(header)

# Install procedure for the header files
#----------------------------------------------------

# ----------------------------------------------------
# --- LIBS -------------------------------------------
# ----------------------------------------------------

#INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include)

ADD_LIBRARY(${LIBRARY_NAME}
  SHARED
  ${${PROJECT_NAME}_SRC})
TARGET_LINK_LIBRARIES(${LIBRARY_NAME} ${Boost_SIGNALS_LIBRARY})

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


