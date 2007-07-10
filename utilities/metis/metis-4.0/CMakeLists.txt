PROJECT(metis_proj)

SET(ARCH_OS "Linux" CACHE STRING "Operating system: either Linux or Linux64")

OPTION(BUILD_SHARED_LIBS "Build shared libraries" OFF)

SUBDIRS(Lib Programs Test)

INSTALL_FILES(/include FILES metis.hpp)
