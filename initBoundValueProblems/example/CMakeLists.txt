##
#  CMake script
##

# Set the name of the project and target:
SET(TARGET "main")
SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

SET(ENV{DEAL_II_DIR} "/nfs/mcfs_comp/home/rudraa/software/dealii-8.4.1/installWithTrilinos")
FIND_PACKAGE(deal.II 8.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

#DEAL_II_INITIALIZE_CACHED_VARIABLES()
#PROJECT(${TARGET} CXX)
#DEAL_II_INVOKE_AUTOPILOT()




FILE(GLOB IGAbase_SRC "../../src/IGA/IGAbase/*.cc" "../../src/IGA/solver/superLU_solver.cc"  "../../src/supplementary/*.cc")
ADD_LIBRARY(IGAbase SHARED ${IGAbase_SRC})
DEAL_II_SETUP_TARGET(IGAbase)

FILE(GLOB model_SRC "../../src/model/*.cc")
ADD_LIBRARY(model SHARED ${model_SRC})
DEAL_II_SETUP_TARGET(model)

FILE (GLOB TARGET_SRC "../../src/IGA/*.cc" "../../src/IGA/function/*.cc" "main.cc")

#FILE(GLOB lib_SRC "../../src/model/*.cc" "../../src/IGA/*.cc""../../src/IGA/IGAbase/*.cc" "../../src/IGA/solver/superLU_solver.cc"  "../../src/supplementary/*.cc" "../../src/IGA/function/*.cc" )
#ADD_LIBRARY(mysrc SHARED ${lib_SRC})
#DEAL_II_SETUP_TARGET(mysrc)


DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()

#add external libraries

SET (hdfLibrary "/sw/arcts/centos7/hdf5/1.8.16-intel-16.0.3/lib/libhdf5_cpp.so")
SET (hdfLibrary2 "/sw/arcts/centos7/hdf5/1.8.16-intel-16.0.3/lib/libhdf5.a")
SET (hdfInclude "/sw/arcts/centos7/hdf5/1.8.16-intel-16.0.3/include")

SET(szip_LIB "/sw/arcts/centos7/szip/2.1/lib/libsz.so")
SET(szip_INC "/sw/arcts/centos7/szip/2.1/include")

SET (SUPERLU_INC "/home/wzhenlin/software/SuperLU_4.3/SRC/")
SET (SUPERLU_LIB "/home/wzhenlin/software/SuperLU_4.3/lib/libsuperlu_4.3.a")

INCLUDE_DIRECTORIES(${hdfInclude} ${SUPERLU_INC} ${szip_INC})
TARGET_LINK_LIBRARIES(main IGAbase model ${hdfLibrary} ${hdfLibrary2} ${SUPERLU_LIB} ${szip_LIB})








#SET (hdfLibrary "/sw/arcts/centos7/hdf5/1.8.16-intel-16.0.3/lib/libhdf5_cpp.so")
#SET (hdfInclude "/sw/arcts/centos7/hdf5/1.8.16-intel-16.0.3/include")
#SET (hdfLibraryDIR "/sw/arcts/centos7/hdf5/1.8.16-intel-16.0.3/lib")

#SET (SUPERLU_INC "/home/wzhenlin/software/SuperLU_4.3/SRC/")
#SET (SUPERLU_LIB "/home/wzhenlin/software/SuperLU_4.3/lib/libsuperlu_4.3.a")
#SET (SUPERLU_LIBDIR "/home/wzhenlin/software/SuperLU_4.3/lib/")
#INCLUDE_DIRECTORIES(${hdfInclude} ${SUPERLU_INC})
#LINK_DIRECTORIES(${hdfLibraryDIR} ${SUPERLU_LIBDIR})
#TARGET_LINK_LIBRARIES(${TARGET} IGAbase model ${hdfLibrary} ${SUPERLU_LIB})


SET (hdfLibrary "/sw/arcts/centos7/hdf5/1.8.16-intel-16.0.3/lib/libhdf5_cpp.so")
SET (hdfInclude "/sw/arcts/centos7/hdf5/1.8.16-intel-16.0.3/include")

SET (SUPERLU_INC "/home/wzhenlin/software/SuperLU_4.3/SRC/")
SET (SUPERLU_LIB "/home/wzhenlin/software/SuperLU_4.3/lib/libsuperlu_4.3.a")
INCLUDE_DIRECTORIES(${hdfInclude} ${SUPERLU_INC})

ADD_LIBRARY(hdf5 STATIC ${hdfLibrary})
ADD_LIBRARY(superLU STATIC ${SUPERLU_LIB})
TARGET_LINK_LIBRARIES(${TARGET} IGAbase model hdf5 superLU)






