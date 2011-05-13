#
#  Copyright
#


ADD_DEFINITIONS(-DDEBUG=2)
IF(WIN32)
	ADD_DEFINITIONS("-DNOMINMAX")
ENDIF(WIN32)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include)
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

# --- tsubmatrix ---
SET(EXECUTABLE_NAME tsubmatrix)
ADD_EXECUTABLE(${EXECUTABLE_NAME}
  tsubmatrix.cpp)
LIST (APPEND executable_list ${EXECUTABLE_NAME})


SET(EXECUTABLE_NAME tsolve)
ADD_EXECUTABLE(${EXECUTABLE_NAME}
  tsolve.cpp)
LIST (APPEND executable_list ${EXECUTABLE_NAME})

# --- tbasic ---
SET(EXECUTABLE_NAME tbasic)
ADD_EXECUTABLE(${EXECUTABLE_NAME}
  tbasic.cpp)
LIST (APPEND executable_list ${EXECUTABLE_NAME})

# --- tstage ---
SET(EXECUTABLE_NAME tstage)
ADD_EXECUTABLE(${EXECUTABLE_NAME}
  tstage.cpp RandomGenerator.cpp MatrixRnd.cpp)
LIST (APPEND executable_list ${EXECUTABLE_NAME})

# --- tdestructiveQR ---
SET(EXECUTABLE_NAME tdestructiveQR)
ADD_EXECUTABLE(${EXECUTABLE_NAME}
  tdestructiveQR.cpp)
LIST (APPEND executable_list ${EXECUTABLE_NAME})

# --- thcod ---
SET(EXECUTABLE_NAME thcod)
ADD_EXECUTABLE(${EXECUTABLE_NAME}
  thcod.cpp MatrixRnd.cpp RandomGenerator.cpp)
LIST (APPEND executable_list ${EXECUTABLE_NAME})

# --- texhcod ---
SET(EXECUTABLE_NAME texhcod)
ADD_EXECUTABLE(${EXECUTABLE_NAME}
  texhcod.cpp MatrixRnd.cpp  gettimeofday.cpp)
LIST (APPEND executable_list ${EXECUTABLE_NAME})

# --- tgivens ---
SET(EXECUTABLE_NAME tgivens)
ADD_EXECUTABLE(${EXECUTABLE_NAME}
  tgivens.cpp)
LIST (APPEND executable_list ${EXECUTABLE_NAME})

# --- tasearch ---
SET(EXECUTABLE_NAME tasearch)
ADD_EXECUTABLE(${EXECUTABLE_NAME}
  tasearch.cpp MatrixRnd.cpp RandomGenerator.cpp)
LIST (APPEND executable_list ${EXECUTABLE_NAME})

# --- trandom ---
SET(EXECUTABLE_NAME trandom)
ADD_EXECUTABLE(${EXECUTABLE_NAME}
  trandom.cpp MatrixRnd.cpp RandomGenerator.cpp gettimeofday.cpp)
LIST (APPEND executable_list ${EXECUTABLE_NAME})
IF(NOT WIN32)
TARGET_LINK_LIBRARIES(${EXECUTABLE_NAME} ${Boost_PROGRAM_OPTIONS_LIBRARY})
ENDIF(NOT WIN32)

# --- tdamping ---
SET(EXECUTABLE_NAME tdamping)
ADD_EXECUTABLE(${EXECUTABLE_NAME}
  tdamping.cpp MatrixRnd.cpp RandomGenerator.cpp gettimeofday.cpp)
LIST (APPEND executable_list ${EXECUTABLE_NAME})
IF(NOT WIN32)
TARGET_LINK_LIBRARIES(${EXECUTABLE_NAME} ${Boost_PROGRAM_OPTIONS_LIBRARY})
ENDIF(NOT WIN32)

# --- tchrono ---
SET(EXECUTABLE_NAME tchrono)
ADD_EXECUTABLE(${EXECUTABLE_NAME}
  tchrono.cpp MatrixRnd.cpp RandomGenerator.cpp gettimeofday.cpp)
LIST (APPEND executable_list ${EXECUTABLE_NAME})

# --- tcompatibility ---
SET(EXECUTABLE_NAME tcompatibility)
ADD_EXECUTABLE(${EXECUTABLE_NAME}
  tcompatibility.cpp MatrixRnd.cpp RandomGenerator.cpp gettimeofday.cpp)
LIST (APPEND executable_list ${EXECUTABLE_NAME})

# --- ttrack ---
SET(EXECUTABLE_NAME ttrack)
ADD_EXECUTABLE(${EXECUTABLE_NAME}
  ttrack.cpp MatrixRnd.cpp RandomGenerator.cpp gettimeofday.cpp)
LIST (APPEND executable_list ${EXECUTABLE_NAME})

# --- tmanu ---
SET(EXECUTABLE_NAME tmanu)
ADD_EXECUTABLE(${EXECUTABLE_NAME}
  tmanu.cpp)
LIST (APPEND executable_list ${EXECUTABLE_NAME})

foreach (exec ${executable_list})
	ADD_DEPENDENCIES(${exec} ${PROJECT_NAME})
	TARGET_LINK_LIBRARIES(${exec} ${PROJECT_NAME})
	LIST2STRING(_cf "${${PROJECT_NAME}_src_CXXFLAGS}")
	LIST2STRING(_lf "${${PROJECT_NAME}_src_LDFLAGS}")
	SET_TARGET_PROPERTIES(${exec}
	  PROPERTIES
	  COMPILE_FLAGS "${_cf}"
	  LINK_FLAGS "${_lf}"
	)
	
	ADD_TEST(${exec} ${exec})
endforeach(exec)
