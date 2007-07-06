PROJECT(metis_lib)

SET(METIS_LIB_SRC coarsen.c fm.c initpart.c match.c ccgraph.c memory.c 
    pmetis.c pqueue.c refine.c util.c timing.c debug.c 
    bucketsort.c graph.c stat.c kmetis.c kwayrefine.c 
    kwayfm.c balance.c ometis.c srefine.c sfm.c separator.c 
    mincover.c mmd.c mesh.c meshpart.c frename.c fortran.c 
    myqsort.c compress.c parmetis.c estmem.c 
    mpmetis.c mcoarsen.c mmatch.c minitpart.c mbalance.c 
    mrefine.c mutil.c mfm.c mkmetis.c mkwayrefine.c mkwayfmh.c 
    mrefine2.c minitpart2.c mbalance2.c mfm2.c 
    kvmetis.c kwayvolrefine.c kwayvolfm.c subdomains.c
)

SET(METIS_LIB_HDR defs.h  macros.h  metis.h  proto.h  rename.h  struct.h)

ADD_LIBRARY(metis ${METIS_LIB_SRC} ${METIS_LIB_HDR})
TARGET_LINK_LIBRARIES(metis m)

INSTALL_FILES(/eng_sci/cubit/metis/metis-4.0/include FILES ${METIS_LIB_HDR})

INSTALL_TARGETS(/eng_sci/cubit/metis/metis-4.0/lib metis)
