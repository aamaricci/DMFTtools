file(GLOB DMFT_EKIN_LOCAL_SRC
  DMFT_EKIN.f90
  )
SET(DMFT_EKIN_src ${DMFT_EKIN_LOCAL_SRC})
ADD_LIBRARY(DMFT_EKINLIB OBJECT ${DMFT_EKIN_src})
ADD_DEPENDENCIES(DMFT_EKINLIB DMFT_CTRL_VARSLIB)
