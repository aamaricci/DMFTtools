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

# Check that Fortran 90 is supported
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
   MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)


IF( (${CMAKE_Fortran_COMPILER_ID} MATCHES Intel) OR (${CMAKE_Fortran_COMPILER_ID} MATCHES GNU))
  MESSAGE(STATUS "Fortran Compiler id = ${CMAKE_Fortran_COMPILER_ID}")
ELSEIF()
  MESSAGE(FATAL_ERROR "Unsupported Fortran compiler (use Intel or GNU). Try export FC=<your FC compiler> ")
ENDIF()  


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
# reSet compiler to MPI is required 
IF(USE_MPI)
  FIND_PACKAGE(MPI REQUIRED)
  IF(MPI_Fortran_FOUND)
    SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
    MESSAGE(STATUS "${Yellow}Set Fortran compiler FC to ${ColourReset}${CMAKE_Fortran_COMPILER}, ID=${CMAKE_Fortran_COMPILER_ID}")
    SET(MPI_CPP "MPI")		#pre-processing option
  ELSE()
    MESSAGE(FATAL_ERROR "MPI Found but No MPI-Fortran compiler can be determined.")    
  ENDIF()
ELSE(USE_MPI)
  SET(MPI_CPP "")
ENDIF(USE_MPI)






############################################################
# Check for .git and extract the current SHA1 revision.
############################################################
SET(GIT_SHA1 0)
SET(GIT_BRANCH "master")
SET(GIT_TAG ${VERSION})
IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  SET(GIT_REPO TRUE)
  FIND_PACKAGE(Git)
  IF(GIT_FOUND)
    EXECUTE_PROCESS(
      COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE GIT_SHA1
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    MESSAGE( STATUS "${Yellow}Git SHA1: ${ColourReset}${GIT_SHA1}" )
    #
    EXECUTE_PROCESS(
      COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE GIT_BRANCH
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    MESSAGE( STATUS "${Yellow}Git branch: ${ColourReset}${GIT_BRANCH}" )
    #
    EXECUTE_PROCESS(
      COMMAND ${GIT_EXECUTABLE} describe --tags
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE GIT_TAG
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    MESSAGE(STATUS "${Yellow}Git tag: ${ColourReset}${GIT_TAG}")
  ENDIF()
ELSE()
  SET(GIT_REPO FALSE)
ENDIF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
SET(VERSION ${GIT_TAG})    	#Update version





############################################################
# Set Fortran options based on BUILD_TYPE and FC ID
############################################################
# -mcmodel=large  this is to remove the 2Gb limit of virtual memory allocation
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU) # this is gfortran
  SET(CMAKE_Fortran_MODDIR_FLAG   "-J")
  SET(CMAKE_Fortran_FLAGS         "-cpp -ffree-line-length-none -fPIC -w ")
   IF(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 10.0.0)
    SET(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -Wno-argument-mismatch")
  ELSE()
    SET(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
  ENDIF()
  SET(CMAKE_Fortran_FLAGS_TESTING "-O2 -funroll-loops")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -p -g -Wsurprising -Waliasing -fwhole-file -fcheck=all -fbacktrace -fbounds-check")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-O3   -funroll-loops")   
elseif(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  SET(CMAKE_Fortran_MODDIR_FLAG   "-module ") #remember the ending white space here 
  SET(CMAKE_Fortran_FLAGS         "-fpp")
  SET(CMAKE_Fortran_FLAGS_TESTING "-O2 -ftz")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-p -O0 -g -fpe1 -warn -debug extended -traceback -check all,noarg_temp_created")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ftz")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES G95)
  SET(CMAKE_Fortran_MODDIR_FLAG   "-fmod=")
  SET(CMAKE_Fortran_FLAGS         "-cpp")
  SET(CMAKE_Fortran_FLAGS_TESTING "-O1  -fsloppy-char")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -ftrace=full")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  
elseif(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
  SET(CMAKE_Fortran_MODDIR_FLAG   "-module ")
  SET(CMAKE_Fortran_FLAGS         "-")
  SET(CMAKE_Fortran_FLAGS         "")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -Mframe")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-O3 -mcmodel=medium -fast -Munroll")
endif()

IF( "${BUILD_TYPE}" MATCHES "DEBUG")
  MESSAGE(STATUS "Fortran Compiler options = ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
ELSEIF("${BUILD_TYPE}" MATCHES "TESTING")
  MESSAGE(STATUS "Fortran Compiler options = ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_TESTING}")
ELSEIF("${BUILD_TYPE}" MATCHES "RELEASE")
  MESSAGE(STATUS "Fortran Compiler options = ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE}")
ENDIF()


#set MPI preprocessing defition:
IF(USE_MPI)
  ADD_DEFINITIONS(-D_MPI)
ELSE(USE_MPI)
  ADD_DEFINITIONS(-D_)
ENDIF(USE_MPI)

IF( "${BUILD_TYPE}" MATCHES "DEBUG")
  ADD_DEFINITIONS(-D_DEBUG)
ENDIF()


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
# Define the actual files and folders that make up the build
############################################################
SET(USER_HOME $ENV{HOME})
SET(USER $ENV{USER})

STRING(TOLOWER "${CMAKE_Fortran_COMPILER_ID}" FC_ID)
STRING(TOLOWER "${CMAKE_BUILD_TYPE}" BUILD_TYPE)

#default prefix is $HOME/opt/<dmft_tools>/<fc_id>/[<git_branch>[/<debug>]]/<version>
SET(PREFIX_DEF_LOC "$ENV{HOME}/opt")
SET(PREFIX_PROJ "${PROJECT_NAME}")
SET(PREFIX_PATH "${FC_ID}")
#if not master branch, include simplified branch name
IF( (NOT GIT_BRANCH MATCHES "master") )
    SET(PREFIX_PATH  "${PREFIX_PATH}/${GIT_BRANCH}")
ENDIF()

#If DEBUG, add /debug to 
IF("${CMAKE_BUILD_TYPE}" MATCHES "DEBUG")
  SET(PREFIX_PATH  "${PREFIX_PATH}/${BUILD_TYPE}")
ENDIF()

#set default prefix:
SET(FULL_VER "${PREFIX_PATH}/${VERSION}")
SET(VERSION_PATH "${PREFIX_PROJ}/${PREFIX_PATH}")
SET(DEFAULT_PREFIX "${PREFIX_DEF_LOC}/${PREFIX_PROJ}/${PREFIX_PATH}/${VERSION}")
SET(PREFIX  "${DEFAULT_PREFIX}" CACHE PATH "Prefix prepended to install directories")
SET(CMAKE_INSTALL_PREFIX "${PREFIX}" CACHE INTERNAL "Prefix prepended to install directories" FORCE)

#Now set module name corresponding to given prefix:
#if prefix is default nothing to be done.
#if user defined prefix prepend
IF(PREFIX STREQUAL "${DEFAULT_PREFIX}")
  SET(MODULE_NAME "${PREFIX_PROJ}/${PREFIX_PATH}/${VERSION}")
ELSE()
  SET(MODULE_NAME "${PREFIX_PROJ}/user_prefix/${PREFIX_PATH}/${VERSION}")
ENDIF()

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
SET(DT_TMP_VER ${CMAKE_BINARY_DIR}/version)
SET(DT_TARGET_LIB ${CMAKE_INSTALL_PREFIX}/lib)
SET(DT_TARGET_INC ${CMAKE_INSTALL_PREFIX}/include)
SET(DT_TARGET_ETC ${CMAKE_INSTALL_PREFIX}/etc)
SET(DT_TARGET_BIN ${CMAKE_INSTALL_PREFIX}/bin)
SET(DT_TARGET_DIR ${CMAKE_INSTALL_PREFIX})

# Have the .mod files placed in the INC folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${DT_TMP_INC})
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_Fortran_MODULE_DIRECTORY})
INCLUDE_DIRECTORIES(BEFORE ${SCIFOR_INCLUDE_DIRS})



SET(SRC_DIR_DMFT_CTRL_VARS ${DT_SRC}/DMFT_CTRL_VARS)
ADD_SUBDIRECTORY(${SRC_DIR_DMFT_CTRL_VARS})

SET(SRC_DIR_DMFT_GF ${DT_SRC}/DMFT_GF)
ADD_SUBDIRECTORY(${SRC_DIR_DMFT_GF})

SET(SRC_DIR_DMFT_GFIO ${DT_SRC}/DMFT_GFIO)
ADD_SUBDIRECTORY(${SRC_DIR_DMFT_GFIO})

SET(SRC_DIR_DMFT_WEISS_FIELD ${DT_SRC}/DMFT_WEISS_FIELD)
ADD_SUBDIRECTORY(${SRC_DIR_DMFT_WEISS_FIELD})

SET(SRC_DIR_DMFT_EKIN ${DT_SRC}/DMFT_EKIN)
ADD_SUBDIRECTORY(${SRC_DIR_DMFT_EKIN})

SET(SRC_DIR_DMFT_TIGHT_BINDING ${DT_SRC}/DMFT_TIGHT_BINDING)
ADD_SUBDIRECTORY(${SRC_DIR_DMFT_TIGHT_BINDING})

SET(SRC_DIR_DMFT_FFT ${DT_SRC}/DMFT_FFT)
ADD_SUBDIRECTORY(${SRC_DIR_DMFT_FFT})

SET(SRC_DIR_DMFT_CONVERGENCE ${DT_SRC}/DMFT_CONVERGENCE)
ADD_SUBDIRECTORY(${SRC_DIR_DMFT_CONVERGENCE})

SET(SRC_DIR_DMFT_MISC ${DT_SRC}/DMFT_MISC)
ADD_SUBDIRECTORY(${SRC_DIR_DMFT_MISC})


ADD_SUBDIRECTORY(${DT_SRC})

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

SET(DT_VERSION_FILE ${DT_SRC}/${PROJECT_NAME}_version_current_sha1.inc)
ADD_CUSTOM_COMMAND(OUTPUT ${DT_VERSION_FILE}
  COMMAND ${CMAKE_COMMAND} -E echo "character(len=41),parameter,public :: ${PROJECT_NAME}_version_sha1 = \"${GIT_SHA1}\"" > ${DT_VERSION_FILE}
  COMMENT "${Yellow}Creating ${DT_VERSION_FILE}: ${ColourReset}"
  VERBATIM)
ADD_CUSTOM_TARGET(DMFT_VERSION_SHA1 DEPENDS ${DT_VERSION_FILE})


ADD_DEPENDENCIES(DMFT_TOOLSLIB
  DMFT_VERSION_SHA1
  DMFT_CTRL_VARSLIB
  DMFT_GFLIB
  DMFT_GFIOLIB
  DMFT_WEISS_FIELDLIB
  DMFT_EKINLIB
  DMFT_TIGHT_BINDINGLIB
  DMFT_FFTLIB
  DMFT_CONVERGENCELIB
  DMFT_MISCLIB
  DMFT_TOOLSLIB
  )



SET(PERMISSION_777
  OWNER_WRITE OWNER_READ OWNER_EXECUTE
  GROUP_WRITE GROUP_READ GROUP_EXECUTE
  WORLD_WRITE WORLD_READ WORLD_EXECUTE)

#REMOVE the entire directory before re-generating it
FILE(REMOVE_RECURSE ${DT_TMP_ETC}/modules)

#Create default .version file if following conditions are verified:
#1. master branch
#2. no debug
#3. default prefix
IF( (GIT_BRANCH MATCHES "master")
    AND (NOT "${CMAKE_BUILD_TYPE}" MATCHES "DEBUG")
    AND (PREFIX STREQUAL "${DEFAULT_PREFIX}")
    ) 
  SET(TMP_VER_MODULE_FILE ${DT_TMP_ETC}/modules/${VERSION_PATH}/.version)
  CONFIGURE_FILE(${DT_ENV}/version.in ${TMP_VER_MODULE_FILE}  @ONLY)
ENDIF()

#Build the DMFT_TOOLS module for environment modules
SET(TMP_MODULE_NAME "${MODULE_NAME}" CACHE PATH "Prefix prepended to install directories")
SET(TMP_ENV_MODULE_FILE ${DT_TMP_ETC}/modules/${TMP_MODULE_NAME})
CONFIGURE_FILE(${DT_ENV}/module.in ${TMP_ENV_MODULE_FILE} @ONLY)


#Build the PKG-CONFIG file
SET(TMP_PKCONFIG_FILE ${DT_TMP_ETC}/${PROJECT_NAME}.pc)
CONFIGURE_FILE( ${DT_ETC}/${PROJECT_NAME}.pc.in ${TMP_PKCONFIG_FILE} @ONLY)

#Build the user CONFIG scripts (sourced in user shell config file, i.e. .bashrc)
SET(CONFIG_SCRIPT dmft_tools_config_user.sh)
SET(TMP_CONFIGVARS_USER_FILE ${DT_TMP_ETC}/${CONFIG_SCRIPT})
CONFIGURE_FILE(${DT_ETC}/dmft_tools_config_user.sh.in ${TMP_CONFIGVARS_USER_FILE} @ONLY)

SET(CONFIG_SCRIPT dmft_tools_config_global.sh)
SET(TMP_CONFIGVARS_GLOBAL_FILE ${DT_TMP_ETC}/${CONFIG_SCRIPT})
CONFIGURE_FILE(${DT_ETC}/dmft_tools_config_global.sh.in ${TMP_CONFIGVARS_GLOBAL_FILE} @ONLY)


FILE(WRITE  ${DT_TMP_VER}  "${VERSION}\n")
INSTALL(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION ${DT_TARGET_INC})

INSTALL(CODE "execute_process(COMMAND \"${CMAKE_COMMAND}\" -E remove ${DT_VERSION_FILE})")

INSTALL(TARGETS dmft_tools DESTINATION ${DT_TARGET_LIB})

INSTALL(DIRECTORY ${DT_TMP_ETC}/ DESTINATION ${DT_TARGET_ETC})

INSTALL(FILES ${TMP_CONFIGVARS_USER_FILE} DESTINATION ${DT_TARGET_BIN}/
  PERMISSIONS ${PERMISSION_777} SETUID)

INSTALL(FILES ${TMP_CONFIGVARS_GLOBAL_FILE} DESTINATION ${DT_TARGET_BIN}/
  PERMISSIONS ${PERMISSION_777} SETUID)

INSTALL(FILES ${DT_TMP_VER} DESTINATION ${DT_TARGET_DIR} 
  PERMISSIONS ${PERMISSION_777} SETUID)

INSTALL(FILES ${DT_TARGET_ETC}/${PROJECT_NAME}.pc DESTINATION $ENV{HOME}/.pkgconfig.d/
  PERMISSIONS ${PERMISSION_777} SETUID)

INSTALL(DIRECTORY ${DT_TARGET_ETC}/modules/ DESTINATION $ENV{HOME}/.modules.d)

MESSAGE( STATUS "${Red}Library version:${ColourReset} ${VERSION}")
MESSAGE( STATUS "${Red}Library will be installed in:${ColourReset} ${CMAKE_INSTALL_PREFIX}")

MESSAGE( STATUS "
>> ${Red}TO CONCLUDE INSTALLATION${ColourReset} <<
Compile with:
$ make

Install with:
$ make install

Uninstall with:
$ make uninstall
")



INSTALL(CODE "MESSAGE(
\"
ADD LIBRARY TO YOUR SYSTEM: 
Pick ONE method below  [add it in your bash profile to make it permanent, e.g. ~/.bashrc]:
${Yellow}Method 1: use the provided ${PROJECT_NAME} environment module${ColourReset}:
   $ module use $HOME/.modules.d
   $ module load ${PROJECT_NAME}/${FC_PLAT}

${Yellow}Method 2: source the config script${ColourReset}:
   $ source ${DT_TARGET_BIN}/dmft_tools_config_user.sh

${Yellow}Method 3: use pkg-config with the provided ${PROJECT_NAME}.pc${ColourReset}:
   $ export PKG_CONFIG_PATH=${DT_TARGET_ETC}/:$PKG_CONFIG_PATH
   $ pkg-config --cflags --libs ${PROJECT_NAME}

${Yellow}Method ADMIN: Add this line to the system shell configuration file, e.g. /etc/bash.bashrc${ColourReset}
   $ source ${DT_TARGET_BIN}dmft_tools_config_global.sh
\")
")


# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_MODULE_PATH}/DistClean.cmake
)



# Uninstall target
if(NOT TARGET uninstall)
  CONFIGURE_FILE(
    "${CMAKE_MODULE_PATH}/ConfigUninstall.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  ADD_CUSTOM_TARGET(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
    "${PROJECT_NAME}/${TMP_MODULE_NAME}" "${PROJECT_NAME}.pc" )
  
ENDIF()


