##################################################
# DMFT_TOOLS PROJECT
##################################################
CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0)
PROJECT(dmft_tools Fortran)
SET(VERSION 2.2.3)

MESSAGE(STATUS "OS: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_VERSION}")
SET(BUILD_TYPE "RELEASE" CACHE STRING "Build type")
STRING(TOUPPER "${BUILD_TYPE}" BUILD_TYPE)
MESSAGE(STATUS "BUILD_TYPE=${BUILD_TYPE}")
SET(CMAKE_BUILD_TYPE ${BUILD_TYPE} CACHE INTERNAL "Build type" FORCE)

# Add our local CMAKE modules to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
INCLUDE(${CMAKE_MODULE_PATH}/ColorsMsg.cmake)


# Set some options the user may choose
OPTION(USE_MPI "Use the MPI library for parallelization" ON)
OPTION(VERBOSE "Build a verbose Makefile" OFF)
IF(VERBOSE)
  SET(CMAKE_VERBOSE_MAKEFILE ON)
ELSE(VERBOSE)
  SET(CMAKE_VERBOSE_MAKEFILE OFF)
ENDIF(VERBOSE)


############################################################
# Look for MPI if needed
############################################################
INCLUDE(${CMAKE_MODULE_PATH}/MpiConfig.cmake)


############################################################
# Find Scifor and add to Include dir
############################################################
FIND_PACKAGE(PkgConfig)
PKG_SEARCH_MODULE(SCIFOR REQUIRED scifor)
IF(${SCIFOR_FOUND})
  MESSAGE(STATUS "${Yellow}Scifor found at: ${ColourReset}${SCIFOR_PREFIX}")
ELSE()
  MESSAGE( FATAL_ERROR "
${Red}Can not find SciFortran in this system. Check SciFor installation. Try adding the provided PKG-CONFIG file (usually in the path-to-scifor/etc) to the environment variable PKG_CONFIG_PATH:${ColourReset}
$ export PKG_CONFIG_PATH=<path_to_scifor_dir>/etc:$PKG_CONFIG_PATH"
    )
ENDIF()



############################################################
# Configure Fortran Compiler Main options
############################################################
INCLUDE(${CMAKE_MODULE_PATH}/FortranCompilerConfig.cmake)




############################################################
# Check for .git and extract the current SHA1 revision.
############################################################
INCLUDE(${CMAKE_MODULE_PATH}/GitConfig.cmake)




############################################################
# Define the actual files and folders that make up the build
############################################################
INCLUDE(${CMAKE_MODULE_PATH}/MainConfig.cmake)




############################################################
############################################################
############################################################
# DEFINE HERE THE DETAILS OF THE LIBRARY CONFIGURATION 


# Define some directories
SET(LIB_SRC ${CMAKE_SOURCE_DIR}/src)
SET(LIB_ETC ${CMAKE_SOURCE_DIR}/etc)
SET(LIB_ENV ${CMAKE_SOURCE_DIR}/etc/environment_modules)
SET(LIB_BIN ${CMAKE_SOURCE_DIR}/bin)
SET(LIB_TMP_INC ${CMAKE_BINARY_DIR}/include)
SET(LIB_TMP_ETC ${CMAKE_BINARY_DIR}/etc)
SET(LIB_TMP_BIN ${CMAKE_BINARY_DIR}/bin)
SET(LIB_TMP_VER ${CMAKE_BINARY_DIR}/version)
SET(LIB_TARGET_LIB ${CMAKE_INSTALL_PREFIX}/lib)
SET(LIB_TARGET_INC ${CMAKE_INSTALL_PREFIX}/include)
SET(LIB_TARGET_ETC ${CMAKE_INSTALL_PREFIX}/etc)
SET(LIB_TARGET_BIN ${CMAKE_INSTALL_PREFIX}/bin)
SET(LIB_TARGET_DIR ${CMAKE_INSTALL_PREFIX})

# Have the .mod files placed in the INC folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${LIB_TMP_INC})
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_Fortran_MODULE_DIRECTORY})
INCLUDE_DIRECTORIES(BEFORE ${SCIFOR_INCLUDE_DIRS})



SET(SRC_DIR_DMFT_CTRL_VARS ${LIB_SRC}/DMFT_CTRL_VARS)
ADD_SUBDIRECTORY(${SRC_DIR_DMFT_CTRL_VARS})

SET(SRC_DIR_DMFT_GF ${LIB_SRC}/DMFT_GF)
ADD_SUBDIRECTORY(${SRC_DIR_DMFT_GF})

SET(SRC_DIR_DMFT_GFIO ${LIB_SRC}/DMFT_GFIO)
ADD_SUBDIRECTORY(${SRC_DIR_DMFT_GFIO})

SET(SRC_DIR_DMFT_WEISS_FIELD ${LIB_SRC}/DMFT_WEISS_FIELD)
ADD_SUBDIRECTORY(${SRC_DIR_DMFT_WEISS_FIELD})

SET(SRC_DIR_DMFT_EKIN ${LIB_SRC}/DMFT_EKIN)
ADD_SUBDIRECTORY(${SRC_DIR_DMFT_EKIN})

SET(SRC_DIR_DMFT_TIGHT_BINDING ${LIB_SRC}/DMFT_TIGHT_BINDING)
ADD_SUBDIRECTORY(${SRC_DIR_DMFT_TIGHT_BINDING})

SET(SRC_DIR_DMFT_FFT ${LIB_SRC}/DMFT_FFT)
ADD_SUBDIRECTORY(${SRC_DIR_DMFT_FFT})

SET(SRC_DIR_DMFT_CONVERGENCE ${LIB_SRC}/DMFT_CONVERGENCE)
ADD_SUBDIRECTORY(${SRC_DIR_DMFT_CONVERGENCE})

SET(SRC_DIR_DMFT_MISC ${LIB_SRC}/DMFT_MISC)
ADD_SUBDIRECTORY(${SRC_DIR_DMFT_MISC})


ADD_SUBDIRECTORY(${LIB_SRC})

# # Build the scifor library and have it placed in the target lib folder
ADD_LIBRARY(dmft_tools
  $<TARGET_OBJECTS:DMFT_CTRL_VARSLIB>
  $<TARGET_OBJECTS:DMFT_GFLIB>
  $<TARGET_OBJECTS:DMFT_GFIOLIB>
  $<TARGET_OBJECTS:DMFT_WEISS_FIELDLIB>
  $<TARGET_OBJECTS:DMFT_EKINLIB>
  $<TARGET_OBJECTS:DMFT_TIGHT_BINDINGLIB>
  $<TARGET_OBJECTS:DMFT_FFTLIB>
  $<TARGET_OBJECTS:DMFT_CONVERGENCELIB>
  $<TARGET_OBJECTS:DMFT_MISCLIB>
  $<TARGET_OBJECTS:DMFT_TOOLSLIB>  
  )



ADD_DEPENDENCIES(dmft_tools
  DMFT_CTRL_VARSLIB
  DMFT_GFLIB
  DMFT_GFIOLIB
  DMFT_WEISS_FIELDLIB
  DMFT_EKINLIB
  DMFT_TIGHT_BINDINGLIB
  DMFT_FFTLIB
  DMFT_CONVERGENCELIB
  DMFT_MISCLIB
  )



SET(PERMISSION_777
  OWNER_WRITE OWNER_READ OWNER_EXECUTE
  GROUP_WRITE GROUP_READ GROUP_EXECUTE
  WORLD_WRITE WORLD_READ WORLD_EXECUTE)


#Build the PKG-CONFIG file
SET(TMP_PKCONFIG_FILE ${LIB_TMP_ETC}/${PROJECT_NAME}.pc)
CONFIGURE_FILE( ${LIB_ETC}/${PROJECT_NAME}.pc.in ${TMP_PKCONFIG_FILE} @ONLY)

############################################################
############################################################
############################################################

INCLUDE(${CMAKE_MODULE_PATH}/PostConfig.cmake)
