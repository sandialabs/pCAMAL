PROJECT (pcamal_proto)

SET (PCAMAL_DIR CACHE FILEPATH "PCamal top-level directory.")

SET (PCAMAL_MPI_LIBRARY CACHE FILEPATH
     "MPI library.")

SET (PCAMAL_MPI_CXX_LIBRARY CACHE FILEPATH
     "MPI CXX library.")

SET (PCAMAL_MPI_INCLUDE CACHE FILEPATH
     "MPI include directory.")

SET (CAMAL_DIR ${PCAMAL_DIR}/utilities/camal)
SET (EXODUS_DIR ${PCAMAL_DIR}/utilities/exodus/exodus-4.01)
SET (NETCDF_DIR ${PCAMAL_DIR}/utilities/netcdf-3.4.snl3)

INCLUDE_DIRECTORIES (
	${CAMAL_DIR}/include
	${EXODUS_DIR}/include
	${NETCDF_DIR}/include
	${PCAMAL_MPI_INCLUDE}
	)

LINK_DIRECTORIES (
	${CAMAL_DIR}/lib/Linux64
	${EXODUS_DIR}/lib/Linux64
	${NETCDF_DIR}/lib/Linux64
	)

SET (PCAMAL_SRCS pcamal_proto.cpp
		 PCExodusFile.cpp
		 PCSweepVolume.cpp)

ADD_EXECUTABLE (pcamal_proto ${PCAMAL_SRCS})

TARGET_LINK_LIBRARIES (pcamal_proto 
		      camal_sweep camal_shr camal_util
		      exoIIv2c401 netcdf dl
		      ${PCAMAL_MPI_LIBRARY}
		      ${PCAMAL_MPI_CXX_LIBRARY})