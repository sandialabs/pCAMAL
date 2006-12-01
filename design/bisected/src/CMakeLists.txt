PROJECT(bisected-file)

SET(PCAMAL_DIR /home/mbsteph/pcamal)
SET(EXODUS_DIR ${PCAMAL_DIR}/utilities/exodus/exodus-4.01)
SET(NETCDF_DIR ${PCAMAL_DIR}/utilities/netcdf-3.4.snl3)

SET(SOURCE bisected-file.cpp)

INCLUDE_DIRECTORIES(${EXODUS_DIR}/include ${NETCDF_DIR}/include)
LINK_DIRECTORIES(${EXODUS_DIR}/lib ${NETCDF_DIR}/lib)

ADD_EXECUTABLE(bisected-file ${SOURCE})

TARGET_LINK_LIBRARIES(bisected-file exoIIv2c401 netcdf)
