FUNCTION(BUILD_ENV_MODULE FILE)
  FILE(WRITE  ${FILE}  "#%Module\n")
  FILE(APPEND ${FILE} "set project ${PROJECT_NAME}\n")
  FILE(APPEND ${FILE} "set root ${PREFIX}\n")
  FILE(APPEND ${FILE} "set plat ${FC_ID}\n")
  FILE(APPEND ${FILE} "set version \"${FULL_VER}\"\n")
  FILE(APPEND ${FILE} "set compiler ${CMAKE_Fortran_COMPILER}\n")
  FILE(READ   ${DT_ENV}/module CONTENTS)
  FILE(APPEND ${FILE} "${CONTENTS}")
ENDFUNCTION()

FUNCTION(BUILD_VERSION_MODULE FILE)
  FILE(WRITE  ${FILE}  "#%Module\n")
  FILE(APPEND ${FILE} "module-version \"./${VERSION}\" default \n")
ENDFUNCTION()


FUNCTION(BUILD_CONFIGVARS_USER FILE)
  FILE(READ   ${DT_BIN}/dmft_tools_config_user.sh CONTENTS)
  FILE(WRITE  ${FILE} "${CONTENTS}")
  FILE(APPEND ${FILE} "add_library_to_system ${CMAKE_INSTALL_PREFIX}\n")
ENDFUNCTION()

FUNCTION(BUILD_CONFIGVARS_GLOBAL FILE)
  FILE(READ   ${DT_BIN}/dmft_tools_config_global.sh CONTENTS)
  FILE(WRITE  ${FILE} "${CONTENTS}")
  FILE(APPEND ${FILE} "add_library_to_system ${CMAKE_INSTALL_PREFIX}\n")
ENDFUNCTION()


FUNCTION(BUILD_PKCONFIG FILE)
  FILE(WRITE ${FILE} "prefix=${CMAKE_INSTALL_PREFIX}\n")
  FILE(READ   ${DT_ETC}/${PROJECT_NAME}.pc CONTENTS)
  FILE(APPEND ${FILE} "${CONTENTS}")
  FILE(APPEND ${FILE} "Version:${VERSION}\n")
ENDFUNCTION()
