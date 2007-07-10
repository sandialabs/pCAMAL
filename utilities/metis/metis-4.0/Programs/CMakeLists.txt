PROJECT(programs)

SET(PMETIS_SRC pmetis.c io.c)
SET(KMETIS_SRC kmetis.c io.c)
SET(OEMETIS_SRC oemetis.c io.c smbfactor.c)
SET(ONMETIS_SRC onmetis.c io.c smbfactor.c)
SET(MESH2DUAL_SRC mesh2dual.c io.c)
SET(MESH2NODAL_SRC mesh2nodal.c io.c)
SET(PARTDMESH_SRC partdmesh.c io.c)
SET(PARTNMESH_SRC partnmesh.c io.c)
SET(GRAPHCHK_SRC  graphchk.c io.c)

INCLUDE_DIRECTORIES(${metis_proj_SOURCE_DIR}/Lib)

ADD_EXECUTABLE(pmetis ${PMETIS_SRC})
TARGET_LINK_LIBRARIES(pmetis metis)

ADD_EXECUTABLE(kmetis ${KMETIS_SRC})
TARGET_LINK_LIBRARIES(kmetis metis)

ADD_EXECUTABLE(oemetis ${OEMETIS_SRC})
TARGET_LINK_LIBRARIES(oemetis metis)

ADD_EXECUTABLE(onmetis ${ONMETIS_SRC})
TARGET_LINK_LIBRARIES(onmetis metis)

ADD_EXECUTABLE(mesh2dual ${MESH2DUAL_SRC})
TARGET_LINK_LIBRARIES(mesh2dual metis)

ADD_EXECUTABLE(mesh2nodal ${MESH2NODAL_SRC})
TARGET_LINK_LIBRARIES(mesh2nodal metis)

ADD_EXECUTABLE(partdmesh ${PARTDMESH_SRC})
TARGET_LINK_LIBRARIES(partdmesh metis)

ADD_EXECUTABLE(partnmesh ${PARTNMESH_SRC})
TARGET_LINK_LIBRARIES(partnmesh metis)

ADD_EXECUTABLE(graphchk ${GRAPHCHK_SRC})
TARGET_LINK_LIBRARIES(graphchk metis)

INSTALL_TARGETS(/bin/${ARCH_OS}
                pmetis kmetis oemetis onmetis mesh2dual mesh2nodal 
		partdmesh partnmesh graphchk)
