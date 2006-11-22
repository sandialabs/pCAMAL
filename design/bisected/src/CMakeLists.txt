PROJECT(bisected-file)

SET(CUBIT_DIR /usr/local/eng_sci/cubit)
SET(EXODUS_VER 4.01)
SET(EXODUS_DIR ${CUBIT_DIR}/exodus/exodus-${EXODUS_VER})
SET(NETCDF_DIR ${CUBIT_DIR}/netcdf/netcdf-3.4.snl3)

SET(SOURCE bisected-file.cpp)

INCLUDE_DIRECTORIES(${EXODUS_DIR}/include ${NETCDF_DIR}/include)
LINK_DIRECTORIES(${EXODUS_DIR}/lib ${NETCDF_DIR}/lib)

ADD_EXECUTABLE(bisected-file ${SOURCE})

TARGET_LINK_LIBRARIES(bisected-file exoIIv2c401 netcdf)
