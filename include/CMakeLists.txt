#
#  Copyright 2010 CNRS
#

SET(${PROJECT_NAME}_HEADERS
     soth/ActiveSet.hpp
     soth/Algebra.hpp
     soth/Allocator.hpp
     soth/BaseY.hpp
     soth/BasicStage.hpp
     soth/Bound.hpp
     soth/debug.hpp
     soth/DestructiveColPivQR.hpp
     soth/Givens.hpp
     soth/HCOD.hpp
     soth/solvers.hpp
     soth/Stage.hpp
     soth/SubMatrix.hpp
)

# Recreate correct path for the headers
#--------------------------------------

SET(fullpath_${PROJECT_NAME}_HEADERS)
FOREACH(lHeader ${${PROJECT_NAME}_HEADERS})
  SET(fullpath_${PROJECT_NAME}_HEADERS
    ${fullpath_${PROJECT_NAME}_HEADERS}
    ${lHeader}
  )
ENDFOREACH(lHeader)

# Install procedure for the header files
#----------------------------------------------------
INSTALL(FILES ${fullpath_${PROJECT_NAME}_HEADERS}
	DESTINATION ${CMAKE_INSTALL_PREFIX}/include/${PROJECT_NAME}
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ OWNER_WRITE
)
