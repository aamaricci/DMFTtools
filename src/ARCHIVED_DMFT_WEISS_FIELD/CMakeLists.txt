file(GLOB ARCHIVED_DMFT_WEISS_FIELD_LOCAL_SRC
  SC_COMMON.f90
  SC_WEISS.f90
  SC_DELTA.f90
  SC_GLOBAL.f90
  DMFT_WEISS_FIELD.f90
  )
SET(ARCHIVED_DMFT_WEISS_FIELD_src ${ARCHIVED_DMFT_WEISS_FIELD_LOCAL_SRC})
ADD_LIBRARY(ARCHIVED_DMFT_WEISS_FIELDLIB OBJECT ${ARCHIVED_DMFT_WEISS_FIELD_src})
ADD_DEPENDENCIES(ARCHIVED_DMFT_WEISS_FIELDLIB DMFT_CTRL_VARSLIB)
