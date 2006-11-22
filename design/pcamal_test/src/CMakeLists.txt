PROJECT(pcamal_test)

SET(CAMAL_DIR /home/mbsteph/Cubit/beta/src/camal)
SET(CAMAL_LIB_DIR /home/mbsteph/Cubit/beta/debug)
SET(CUBIT_DIR /usr/local/eng_sci/cubit)
SET(EXODUS_VER 4.01)
SET(EXODUS_DIR ${CUBIT_DIR}/exodus/exodus-${EXODUS_VER})
SET(NETCDF_DIR ${CUBIT_DIR}/netcdf/netcdf-3.4.snl3)

SET(SOURCE pcamal_test.cpp PCExodusFile.cpp PCSweepVolume.cpp)

INCLUDE_DIRECTORIES(
	${CAMAL_DIR}/include
	${EXODUS_DIR}/include
	${NETCDF_DIR}/include
	)
LINK_DIRECTORIES(
	${CAMAL_LIB_DIR}
	${EXODUS_DIR}/lib
	${NETCDF_DIR}/lib
	)

ADD_EXECUTABLE(pcamal_test ${SOURCE})

TARGET_LINK_LIBRARIES(pcamal_test camal_sweep camal_shr cubit_util 
		      exoIIv2c401 netcdf dl)
