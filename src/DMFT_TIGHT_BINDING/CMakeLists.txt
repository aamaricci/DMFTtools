file(GLOB DMFT_TIGHT_BINDING_LOCAL_SRC
  TB_COMMON.f90
  TB_IO.f90
  TB_BASIS.f90
  TB_EFERMI.f90
  TB_BUILD.f90
  TB_SOLVE.f90
  TB_FSURF.f90
  TB_WANNIER90.f90
  DMFT_TIGHT_BINDING.f90
  )
SET(DMFT_TIGHT_BINDING_src ${DMFT_TIGHT_BINDING_LOCAL_SRC})
ADD_LIBRARY(DMFT_TIGHT_BINDINGLIB OBJECT ${DMFT_TIGHT_BINDING_src})
ADD_DEPENDENCIES(DMFT_TIGHT_BINDINGLIB DMFT_CTRL_VARSLIB)
ADD_DEPENDENCIES(DMFT_TIGHT_BINDINGLIB DMFT_GFLIB)
