##################################################
# DMFT_TOOLS PROJECT
##################################################
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
PROJECT(dmft_tools Fortran)
SET(VERSION 1.2.7)

MESSAGE(STATUS "OS: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_VERSION}")
SET(BUILD_TYPE "RELEASE" CACHE STRING "Build type")
STRING(TOUPPER "${BUILD_TYPE}" BUILD_TYPE)
MESSAGE(STATUS "BUILD_TYPE=${BUILD_TYPE}")
SET(CMAKE_BUILD_TYPE ${BUILD_TYPE} CACHE INTERNAL "Build type" FORCE)

# Add our local CMAKE modules to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")

INCLUDE(${CMAKE_MODULE_PATH}/ColorsMsg.cmake)

# Check that Fortran 90 is supported
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
   MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)


IF( (${CMAKE_Fortran_COMPILER_ID} MATCHES Intel) OR (${CMAKE_Fortran_COMPILER_ID} MATCHES GNU))
  MESSAGE(STATUS "Fortran Compiler id = ${CMAKE_Fortran_COMPILER_ID}")
ELSEIF()
  MESSAGE(FATAL_ERROR "Unsupported Fortran compiler (use Intel or GNU) ")
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

# # This executes code that sets the compile flags for DEBUG, RELEASE, and TESTING.
# INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake)


############################################################
# Check for .git and extract the current SHA1 revision.
############################################################
IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  FIND_PACKAGE(Git)
  IF(GIT_FOUND)
    EXECUTE_PROCESS(
      COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE GIT_SHA1
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    MESSAGE( STATUS "${Yellow}Git SHA1: ${ColourReset}${GIT_SHA1}" )

    EXECUTE_PROCESS(
      COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE GIT_BRANCH
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    MESSAGE( STATUS "${Yellow}Git branch: ${ColourReset}${GIT_BRANCH}" )
    
    EXECUTE_PROCESS(
      COMMAND ${GIT_EXECUTABLE} describe --tags
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE GIT_TAG
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    MESSAGE(STATUS "${Yellow}Git tag: ${ColourReset}${GIT_TAG}")
  ELSE(GIT_FOUND)
    SET(GIT_SHA1 0)
    SET(GIT_BRANCH "master")
    SET(GIT_TAG 0.0.0)
  ENDIF(GIT_FOUND)
ENDIF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
SET(VERSION ${GIT_TAG})    	#Update version




############################################################
# Set Fortran options based on BUILD_TYPE and FC ID
############################################################
# -mcmodel=large  this is to remove the 2Gb limit of virtual memory allocation
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU) # this is gfortran
  SET(CMAKE_Fortran_MODDIR_FLAG   "-J")
  SET(CMAKE_Fortran_FLAGS         "-cpp -ffree-line-length-none") 
  SET(CMAKE_Fortran_FLAGS_TESTING "-O2 -funroll-loops")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -p -g -Wsurprising -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -fbounds-check")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-O3  -ffast-math -march=native -funroll-loops")   
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

#set MPI preprocessing defition:
IF(USE_MPI)
  ADD_DEFINITIONS(-D_MPI)
ELSE(USE_MPI)
  ADD_DEFINITIONS(-D_)
ENDIF(USE_MPI)


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


#make sure FC_PLAT is lower case
IF(GIT_BRANCH MATCHES "master")
  STRING(TOLOWER "${CMAKE_Fortran_COMPILER_ID}" FC_PLAT)
ELSE()
  SET(FC_PLAT "${GIT_BRANCH}")
ENDIF()

IF("${CMAKE_BUILD_TYPE}" MATCHES "DEBUG")
  SET(PREFIX "$ENV{HOME}/opt/${PROJECT_NAME}/${FC_PLAT}_${CMAKE_BUILD_TYPE}" CACHE PATH "Prefix prepended to install directories")
ELSE()
  SET(PREFIX "$ENV{HOME}/opt/${PROJECT_NAME}/${FC_PLAT}" CACHE PATH "Prefix prepended to install directories")
ENDIF()
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

SET(SRC_DIR_DMFT_GLOC ${DT_SRC}/DMFT_GLOC)
ADD_SUBDIRECTORY(${SRC_DIR_DMFT_GLOC})

SET(SRC_DIR_DMFT_GK ${DT_SRC}/DMFT_GK)
ADD_SUBDIRECTORY(${SRC_DIR_DMFT_GK})

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

SET(SRC_DIR_DMFT_VECTORS ${DT_SRC}/DMFT_VECTORS)
ADD_SUBDIRECTORY(${SRC_DIR_DMFT_VECTORS})

SET(SRC_DIR_DMFT_CONVERGENCE ${DT_SRC}/DMFT_CONVERGENCE)
ADD_SUBDIRECTORY(${SRC_DIR_DMFT_CONVERGENCE})

SET(SRC_DIR_DMFT_MISC ${DT_SRC}/DMFT_MISC)
ADD_SUBDIRECTORY(${SRC_DIR_DMFT_MISC})


ADD_SUBDIRECTORY(${DT_SRC})

# # Build the scifor library and have it placed in the target lib folder
ADD_LIBRARY(dmft_tools
  $<TARGET_OBJECTS:DMFT_CTRL_VARSLIB>
  $<TARGET_OBJECTS:DMFT_GLOCLIB>
  $<TARGET_OBJECTS:DMFT_GKLIB>
  $<TARGET_OBJECTS:DMFT_GFIOLIB>
  $<TARGET_OBJECTS:DMFT_WEISS_FIELDLIB>
  $<TARGET_OBJECTS:DMFT_EKINLIB>
  $<TARGET_OBJECTS:DMFT_TIGHT_BINDINGLIB>
  $<TARGET_OBJECTS:DMFT_FFTLIB>
  $<TARGET_OBJECTS:DMFT_VECTORSLIB>
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
  DMFT_GLOCLIB
  DMFT_GKLIB
  DMFT_GFIOLIB
  DMFT_WEISS_FIELDLIB
  DMFT_EKINLIB
  DMFT_TIGHT_BINDINGLIB
  DMFT_FFTLIB
  DMFT_VECTORSLIB  
  DMFT_CONVERGENCELIB
  DMFT_MISCLIB
  DMFT_TOOLSLIB
  )


# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake
)




SET(PERMISSION_777
  OWNER_WRITE OWNER_READ OWNER_EXECUTE
  GROUP_WRITE GROUP_READ GROUP_EXECUTE
  WORLD_WRITE WORLD_READ WORLD_EXECUTE)


#Contains the functions used in the following:
#BUILD_ENV_MODULE
#BUILD_CONFIGVARS
#BUILD_PKCONFIG
INCLUDE(${CMAKE_MODULE_PATH}/PostBuildDmftTools.cmake)

#Build the dmft_tools module for environment modules
SET(TMP_ENV_MODULE_FILE ${DT_TMP_ETC}/modules/${PROJECT_NAME}/${FC_PLAT})
BUILD_ENV_MODULE(${TMP_ENV_MODULE_FILE})

#Build the user config script (sourced in user shell config file, i.e. .bashrc)
SET(TMP_CONFIGVARS_USER_FILE ${DT_TMP_BIN}/dmft_tools_config_user.sh)
BUILD_CONFIGVARS_USER(${TMP_CONFIGVARS_USER_FILE})

#Build the global config script (sourced in system shell config file, i.e. /etc/bash.bashrc)
SET(TMP_CONFIGVARS_GLOBAL_FILE ${DT_TMP_BIN}/dmft_tools_config_global.sh)
BUILD_CONFIGVARS_GLOBAL(${TMP_CONFIGVARS_GLOBAL_FILE})

SET(TMP_PKCONFIG_FILE ${DT_TMP_ETC}/${PROJECT_NAME}.pc)
BUILD_PKCONFIG(${TMP_PKCONFIG_FILE})

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


MESSAGE( STATUS "${Red}Library version:${ColourReset} ${VERSION}")
MESSAGE( STATUS "${Red}Library will be installed in:${ColourReset} ${CMAKE_INSTALL_PREFIX}")

MESSAGE( STATUS "${Yellow}To conclude installation:${ColourReset} 
$ make
$ make install
$ make post-install
")




ADD_CUSTOM_TARGET(
  post-install
  COMMAND "${CMAKE_COMMAND}" -P ${CMAKE_MODULE_PATH}/HomeInstall.cmake ${CMAKE_SOURCE_DIR} ${CMAKE_INSTALL_PREFIX} ${PROJECT_NAME} ${FC_PLAT}
  WORKING_DIRECTORY
  "${CMAKE_BINARY_DIR}"
)


