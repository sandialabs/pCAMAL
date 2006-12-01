PROJECT(pcamal_test)

SET(PCAMAL_DIR /home/mbsteph/pcamal)
SET(CAMAL_DIR ${PCAMAL_DIR}/utilities/camal)
SET(EXODUS_DIR ${PCAMAL_DIR}/utilities/exodus/exodus-4.01)
SET(NETCDF_DIR ${PCAMAL_DIR}/utilities/netcdf-3.4.snl3)

SET(SOURCE pcamal_test.cpp PCExodusFile.cpp PCSweepVolume.cpp)

INCLUDE_DIRECTORIES(
	${CAMAL_DIR}/include
	${EXODUS_DIR}/include
	${NETCDF_DIR}/include
	)
LINK_DIRECTORIES(
	${CAMAL_DIR}/lib
	${EXODUS_DIR}/lib
	${NETCDF_DIR}/lib
	)

ADD_EXECUTABLE(pcamal_test ${SOURCE})

TARGET_LINK_LIBRARIES(pcamal_test camal_sweep camal_shr camal_util 
		      exoIIv2c401 netcdf dl)
