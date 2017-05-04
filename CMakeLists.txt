##################################################
# SCIFOR PROJECT (TEST CMAKE(
##################################################
CMAKE_MINIMUM_REQUIRED(VERSION 3.0.2)
PROJECT(dmft_tools Fortran)
SET(VERSION 0.7.0)
SET(CMAKE_BUILD_TYPE TESTING)

# Add our local CMAKE modules to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")

# Check that Fortran 90 is supported
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
   MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

# Set some options the user may choose
OPTION(USE_MPI "Use the MPI library for parallelization" ON)
OPTION(VERBOSE "Build a verbose Makefile" OFF)

# reSet compiler to MPI is required 
IF(USE_MPI)
  FIND_PACKAGE(MPI REQUIRED)
  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
  SET(MPI_CPP "MPI")
  MESSAGE("-+ Set Fortran compiler FC to ${CMAKE_Fortran_COMPILER}")
ELSE(USE_MPI)
  SET(MPI_CPP "")
ENDIF(USE_MPI)

# Set verbosity of the Makefile according to option
IF(VERBOSE)
  SET(CMAKE_VERBOSE_MAKEFILE ON)
ELSE(VERBOSE)
  SET(CMAKE_VERBOSE_MAKEFILE OFF)
ENDIF(VERBOSE)

# This executes code that sets the compile flags for DEBUG, RELEASE, and TESTING.  
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake)

#workaround to set the Module directory output flag for the Intel Fortran
#compiler which should not be -J but -module (with a blank space at the end).
#soon enough I will find a more elegant solution.
IF(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  SET(CMAKE_Fortran_MODDIR_FLAG "-module ")
ENDIF()
MESSAGE("CMAKE_Fortran_MODDIR_FLAG:${CMAKE_Fortran_MODDIR_FLAG}")
#set MPI preprocessing defition:
IF(USE_MPI)
  ADD_DEFINITIONS(-D_MPI)
ELSE(USE_MPI)
  ADD_DEFINITIONS(-D_)
ENDIF(USE_MPI)


############################################################
# Define the actual files and folders that make up the build
############################################################
SET(USER_HOME $ENV{HOME})
SET(SF_ROOT $ENV{SFROOT})
IF(NOT EXISTS ${SF_ROOT})
  message( FATAL_ERROR "Can not locate SciFortran. Check environment variable SF_ROOT or specify it to cmake using -DSF_ROOT=path-to-scifor-root" )
ENDIF()

#make sure FC_PLAT is lowecase
STRING(TOLOWER "${CMAKE_Fortran_COMPILER_ID}" FC_PLAT)

SET(PREFIX  "$ENV{HOME}/opt/${PROJECT_NAME}/${VERSION}/${FC_PLAT}" CACHE PATH "Prefix prepended to install directories")

SET(CMAKE_INSTALL_PREFIX "${PREFIX}" CACHE INTERNAL "Prefix prepended to install directories" FORCE)

# Define the library name
SET(DMFT_TOOLSLIB dmft_tools)

# Define some directories
SET(DT_SRC ${CMAKE_SOURCE_DIR}/src)
SET(DT_ETC ${CMAKE_SOURCE_DIR}/etc)
SET(DT_ENV ${CMAKE_SOURCE_DIR}/etc/environment_modules)
SET(DT_BIN ${CMAKE_SOURCE_DIR}/bin)
SET(DT_TMP_INC ${CMAKE_BINARY_DIR}/include)
SET(DT_TMP_ETC ${CMAKE_BINARY_DIR}/etc)
SET(DT_TMP_BIN ${CMAKE_BINARY_DIR}/bin)
SET(DT_TARGET_LIB ${CMAKE_INSTALL_PREFIX}/lib)
SET(DT_TARGET_INC ${CMAKE_INSTALL_PREFIX}/include)
SET(DT_TARGET_ETC ${CMAKE_INSTALL_PREFIX}/etc)
SET(DT_TARGET_BIN ${CMAKE_INSTALL_PREFIX}/bin)

# Have the .mod files placed in the INC folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${DT_TMP_INC})
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_Fortran_MODULE_DIRECTORY})
INCLUDE_DIRECTORIES(BEFORE ${SF_ROOT}/include)

ADD_SUBDIRECTORY(${DT_SRC})

# # Build the scifor library and have it placed in the target lib folder
ADD_LIBRARY(dmft_tools
  $<TARGET_OBJECTS:DMFT_TOOLSLIB>
  )


# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake
)

INCLUDE(${CMAKE_MODULE_PATH}/PostBuildDmftTools.cmake)

SET(DT_TMP_ENV_MODULE_FILE ${DT_TMP_ETC}/modules/${PROJECT_NAME}/${FC_PLAT})
BUILD_DT_ENV_MODULE(${DT_TMP_ENV_MODULE_FILE})
SET(DT_TMP_CONFIGVARS_FILE ${DT_TMP_BIN}/configvars.sh)
BUILD_DT_CONFIGVARS(${DT_TMP_CONFIGVARS_FILE})

INSTALL(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION ${DT_TARGET_INC})

INSTALL(TARGETS dmft_tools DESTINATION ${DT_TARGET_LIB})

INSTALL(DIRECTORY ${DT_TMP_ETC}/ DESTINATION ${DT_TARGET_ETC})

INSTALL(FILES ${DT_TMP_CONFIGVARS_FILE} DESTINATION ${DT_TARGET_BIN}/
  PERMISSIONS
  OWNER_WRITE OWNER_READ OWNER_EXECUTE
  GROUP_READ GROUP_WRITE GROUP_EXECUTE
  WORLD_READ WORLD_WRITE WORLD_EXECUTE)
