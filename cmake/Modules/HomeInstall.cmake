SET(CMAKE_SOURCE_DIR ${CMAKE_ARGV3})
SET(CMAKE_INSTALL_PREFIX ${CMAKE_ARGV4})
SET(PROJECT_NAME ${CMAKE_ARGV5})
SET(FC_PLAT ${CMAKE_ARGV6})

SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")
SET(DT_TARGET_LIB ${CMAKE_INSTALL_PREFIX}/lib)
SET(DT_TARGET_INC ${CMAKE_INSTALL_PREFIX}/include)
SET(DT_TARGET_ETC ${CMAKE_INSTALL_PREFIX}/etc)
SET(DT_TARGET_BIN ${CMAKE_INSTALL_PREFIX}/bin)
SET(PKCONFIG_FILE ${DT_TARGET_ETC}/${PROJECT_NAME}.pc)


INCLUDE(${CMAKE_MODULE_PATH}/ColorsMsg.cmake)

FILE(INSTALL ${DT_TARGET_ETC}/modules/ DESTINATION $ENV{HOME}/.modules.d
  USE_SOURCE_PERMISSIONS)

FILE(INSTALL ${PKCONFIG_FILE} DESTINATION $ENV{HOME}/.pkgconfig.d
  USE_SOURCE_PERMISSIONS)


MESSAGE( STATUS "After installation usage options:")
MESSAGE( STATUS
  "

>>${Red}GLOBAL INSTALLATION${ColourReset}:
${Yellow}(if you are admin of a multi-user machine and you want users to automatically access the library)${ColourReset}
${Green} Add this line to the system shell configuration file (i.e. /etc/bash.bashrc)${ColourReset}:
   $ source ${DT_TARGET_BIN}/dmft_tools_config_global.sh


>>${Red}USER INSTALLATION${ColourReset}:
${Yellow}(if you are regular user of a machine, any of the followings will give you access to the library)${ColourReset}
Pick ONE choice or add in your bash profile (i.e. ~/.bashrc):
${Green}A. source the config script: ${DT_TARGET_BIN}/dmft_tools_config_user.sh, i.e.${ColourReset}:
   $ source ${DT_TARGET_BIN}/dmft_tools_config_user.sh

${Green}B. use the provided ${PROJECT_NAME} environment module in ${DT_TARGET_ETC}, i.e.${ColourReset}:
   $ module use $HOME/.modules.d
   $ module load ${PROJECT_NAME}/${FC_PLAT}

${Green}C. use pkg-config with the provided ${PROJECT_NAME}.pc in ${DT_TARGET_ETC}, i.e.${ColourReset}:
   $ export PKG_CONFIG_PATH=${DT_TARGET_ETC}/:$PKG_CONFIG_PATH
   $ pkg-config --cflags --libs ${PROJECT_NAME} (to get lib info)
")
