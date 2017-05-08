##################################################
# SCIFOR PROJECT (TEST CMAKE(
##################################################
CMAKE_MINIMUM_REQUIRED(VERSION 3.0.2)
PROJECT(dmft_tools Fortran)
SET(VERSION 0.7.1)
SET(CMAKE_BUILD_TYPE RELEASE)

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
  MESSAGE(STATUS "Set Fortran compiler FC to ${CMAKE_Fortran_COMPILER}")
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
INCLUDE(${CMAKE_MODULE_PATH}/ColorsMsg.cmake)
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
    MESSAGE( STATUS "${Magenta}Git version: ${GIT_SHA1}${ColourReset}" )
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




SET(PERMISSION_777
  OWNER_WRITE OWNER_READ OWNER_EXECUTE
  GROUP_WRITE GROUP_READ GROUP_EXECUTE
  WORLD_WRITE WORLD_READ WORLD_EXECUTE)

SET(PERMISSION_733
  OWNER_WRITE OWNER_READ OWNER_EXECUTE
  GROUP_READ GROUP_EXECUTE
  WORLD_READ WORLD_EXECUTE)

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

# INSTALL(DIRECTORY  ${DT_TMP_ETC}/modules/ DESTINATION $ENV{HOME}/.modules.d/libraries/
#   USE_SOURCE_PERMISSIONS
#   FILE_PERMISSIONS ${PERMISSION_777} SETUID)

# INSTALL(FILES  ${TMP_PKCONFIG_FILE} DESTINATION $ENV{HOME}/.pkgconfig.d/
#   PERMISSIONS ${PERMISSION_733} SETUID)


MESSAGE( STATUS "After installation usage options:")
MESSAGE( STATUS
  "Pick one choice (put it in your bash profile file to make it permanent):
${Magenta}A. source the file: ${SF_TARGET_BIN}/configvars.sh ${ColourReset}
i.e.
$ source ${SF_TARGET_BIN}/configvars.sh

${Magenta}B. use the provided ${PROJECT_NAME} environment module in ${DT_TARGET_ETC}:${ColourReset}
i.e.
$ mkdir -p $HOME/.modules.d
$ module use $HOME/.modules.d
$ cp -r ${DT_TARGET_ETC}/modules $HOME/.modules.d
$ module load ${PROJECT_NAME}/${FC_PLAT}

${Magenta}C. use pkg-config with the provided ${PROJECT_NAME}.pc in ${DT_TARGET_ETC}: ${ColourReset}
i.e.
$ export PKG_CONFIG_PATH=${DT_TARGET_ETC}/:$PKG_CONFIG_PATH
$ pkg-config --cflags --libs ${PROJECT_NAME} (to get lib info)
")
