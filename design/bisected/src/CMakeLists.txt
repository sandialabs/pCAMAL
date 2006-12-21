PROJECT(bisected-file)

SET(PCAMAL_OS CACHE STRING "Operating system (Linux, Linux64)")
SET(PCAMAL_DIR CACHE FILEPATH "pCAMAL top-level directory.")

SET(EXODUS_DIR ${PCAMAL_DIR}/utilities/exodus/exodus-4.01)
SET(NETCDF_DIR ${PCAMAL_DIR}/utilities/netcdf-3.4.snl3)

SET(SOURCE bisected-file.cpp)

INCLUDE_DIRECTORIES(${EXODUS_DIR}/include ${NETCDF_DIR}/include)
LINK_DIRECTORIES(${EXODUS_DIR}/lib/${PCAMAL_OS}
	 ${NETCDF_DIR}/lib/${PCAMAL_OS})

ADD_EXECUTABLE(bisected-file ${SOURCE})

TARGET_LINK_LIBRARIES(bisected-file exoIIv2c401 netcdf)
