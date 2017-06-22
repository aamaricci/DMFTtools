##################################################
# DMFT_TOOLS PROJECT
##################################################
CMAKE_MINIMUM_REQUIRED(VERSION 3.0.2)
PROJECT(dmft_tools Fortran)
SET(VERSION 0.7.2d)

SET(BUILD_TYPE "RELEASE" CACHE STRING "Build type")
MESSAGE(STATUS "BUILD_TYPE=${BUILD_TYPE}")
SET(CMAKE_BUILD_TYPE ${BUILD_TYPE} CACHE INTERNAL "Build type" FORCE)

# Add our local CMAKE modules to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")

INCLUDE(${CMAKE_MODULE_PATH}/ColorsMsg.cmake)

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
  MESSAGE(STATUS "${Yellow}Set Fortran compiler FC to ${ColourReset}${CMAKE_Fortran_COMPILER}")
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


#Check for .git and extract the current SHA1 revision.
#create the current scifor_version.inc file
IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  FIND_PACKAGE(Git)
  IF(GIT_FOUND)
    EXECUTE_PROCESS(
      COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE GIT_SHA1
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    MESSAGE( STATUS "${Yellow}Git version: ${ColourReset}${GIT_SHA1}" )
  ELSE(GIT_FOUND)
    SET(GIT_SHA1 0)
  ENDIF(GIT_FOUND)
ENDIF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
SET(DT_VERSION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/src/${PROJECT_NAME}_version.inc)
FILE(WRITE  ${DT_VERSION_FILE}
  "character(len=41),parameter,public :: ${PROJECT_NAME}_version_sha1 = \"${GIT_SHA1}\"\n")


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
SET(USER $ENV{USER})

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
INCLUDE_DIRECTORIES(BEFORE ${SCIFOR_INCLUDE_DIRS})

ADD_SUBDIRECTORY(${DT_SRC})

# # Build the scifor library and have it placed in the target lib folder
ADD_LIBRARY(dmft_tools
  $<TARGET_OBJECTS:DMFT_TOOLSLIB>
  )


# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake
)




SET(PERMISSION_777
  OWNER_WRITE OWNER_READ OWNER_EXECUTE
  GROUP_WRITE GROUP_READ GROUP_EXECUTE
  WORLD_WRITE WORLD_READ WORLD_EXECUTE)


INCLUDE(${CMAKE_MODULE_PATH}/PostBuildDmftTools.cmake)

SET(TMP_ENV_MODULE_FILE ${DT_TMP_ETC}/modules/${PROJECT_NAME}/${FC_PLAT})
BUILD_ENV_MODULE(${TMP_ENV_MODULE_FILE})

SET(TMP_CONFIGVARS_FILE ${DT_TMP_BIN}/configvars.sh)
BUILD_CONFIGVARS(${TMP_CONFIGVARS_FILE})

SET(TMP_PKCONFIG_FILE ${DT_TMP_ETC}/${PROJECT_NAME}.pc)
BUILD_PKCONFIG(${TMP_PKCONFIG_FILE})


INSTALL(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION ${DT_TARGET_INC})


INSTALL(TARGETS dmft_tools DESTINATION ${DT_TARGET_LIB})

INSTALL(DIRECTORY ${DT_TMP_ETC}/ DESTINATION ${DT_TARGET_ETC})

INSTALL(FILES ${TMP_CONFIGVARS_FILE} DESTINATION ${DT_TARGET_BIN}/
  PERMISSIONS ${PERMISSION_777} SETUID)

MESSAGE( STATUS "${Yellow}To conclude installation:${ColourReset} 
$ make
$ make install
$ make post-install
")


ADD_CUSTOM_TARGET(
  post-install
  COMMAND "${CMAKE_COMMAND}" -P ${CMAKE_MODULE_PATH}/HomeInstall.cmake ${CMAKE_SOURCE_DIR} ${CMAKE_INSTALL_PREFIX} ${PROJECT_NAME}
  WORKING_DIRECTORY
  "${CMAKE_BINARY_DIR}"
)


