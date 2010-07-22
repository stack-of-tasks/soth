#
#  Copyright 2010 CNRS
#

SET(${PROJECT_NAME}_HEADERS
    soth/soth.h
    soth/Algebra.h
    soth/SubMatrix.h
    soth/DestructiveColPivQR.h
    soth/debug.h
    soth/ActiveSet.h
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
	DESTINATION ${CMAKE_INSTALL_PREFIX}/include/dynamic-graph/tutorial
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ OWNER_WRITE
)
