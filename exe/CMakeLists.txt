INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include)
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

# Static library that will avoid useless recompilation.
#add_library(utils STATIC
#  ../unitTesting/gettimeofday.cpp ../unitTesting/RandomGenerator.cpp
#)

# --- tsubmatrix ---
SET(EXECUTABLE_NAME conic_simplification)
ADD_EXECUTABLE(${EXECUTABLE_NAME} 
  conic_simplification.cpp)
ADD_DEPENDENCIES(${EXECUTABLE_NAME} ${PROJECT_NAME})
ADD_DEPENDENCIES(${EXECUTABLE_NAME} utils)
TARGET_LINK_LIBRARIES(${EXECUTABLE_NAME} ${PROJECT_NAME} utils)

# --- tbasic ---
SET(EXECUTABLE_NAME simple)
ADD_EXECUTABLE(${EXECUTABLE_NAME} 
  simple.cpp)
ADD_DEPENDENCIES(${EXECUTABLE_NAME} ${PROJECT_NAME})
ADD_DEPENDENCIES(${EXECUTABLE_NAME} utils)
TARGET_LINK_LIBRARIES(${EXECUTABLE_NAME} ${PROJECT_NAME} utils)

