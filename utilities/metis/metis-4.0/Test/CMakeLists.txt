PROJECT(mtest)

SET(MTEST_SRC mtest.c ${metis_proj_SOURCE_DIR}/Programs/io.c)

INCLUDE_DIRECTORIES(${metis_proj_SOURCE_DIR}/Lib)

ADD_EXECUTABLE(mtest ${MTEST_SRC})
TARGET_LINK_LIBRARIES(mtest metis)
