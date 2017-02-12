##
#  CMake script
##

#Code structure and usage
##---------------------------------------------------------------------------------------------------------------------
#The code is basded on openIGA and contains three level 
#level 1:generic code for initial boundary value problem, you do not need to change anything inside follow folders
#--        defectsElasticity/src  defectsElasticity/include
#
#level 2:generic code for defects modling, usually you do not need to change anything inside follow folder 
#--        defectsElasticity/initBoundValProbs/defects/src
#
#level 3:code for implementing specific defects type, change them as you need 
#--        defectsElasticity/initBoundValProbs/defects/edge_dislication/single/  .....
#
# Usage
# in build folder where you found this CMakeList.txt, change it following guides below
# $cmake CMakeList.txt
# $make run 
##---------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------



# Set the name of the project and target:
SET(TARGET "main")


#openIGA is shared lib
##---------------------------------------------------------------------------------------------------------------------
SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
##---------------------------------------------------------------------------------------------------------------------



#Include deal.II
##---------------------------------------------------------------------------------------------------------------------
#usually DEAL_II_DIR for explicitly deal.II path; deal.II_DIR and DEAL_II_DIR could be preloaded enviromental variable

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
##---------------------------------------------------------------------------------------------------------------------



#Include external lib as necessary (if they cannot be found automatically);
##---------------------------------------------------------------------------------------------------------------------
SET(IGA_INC "/home/snverner/software/openIGA-master/include/")
SET(IGA_LIB "/home/snverner/software/openIGA-master/lib/libopenIGA.so")

SET (hdfLibrary "/sw/arcts/centos7/hdf5/1.8.16-intel-16.0.3/lib/libhdf5_cpp.so")
SET (hdfLibrary2 "/sw/arcts/centos7/hdf5/1.8.16-intel-16.0.3/lib/libhdf5.a")
SET (hdfInclude "/sw/arcts/centos7/hdf5/1.8.16-intel-16.0.3/include")

SET(szip_LIB "/sw/arcts/centos7/szip/2.1/lib/libsz.so")
SET(szip_INC "/sw/arcts/centos7/szip/2.1/include")

SET (SUPERLU_INC "/home/snverner/software/SuperLU_4.3/SRC/")
SET (SUPERLU_LIB "/home/snverner/software/SuperLU_4.3/lib/libsuperlu_4.3.a")
##---------------------------------------------------------------------------------------------------------------------


#Include all source files
##---------------------------------------------------------------------------------------------------------------------
#set shortcut for three directories (two source file folder and one head file folders)
SET(src_DIR "/home/snverner/workspace/Defect2/src")
SET(src_INC "/home/snverner/workspace/Defect2/include")
SET(srcDefect_DIR "/home/snverner/workspace/Defect2/initBoundValProbs/defects/src")


#Include all source files, usually you do NOT need to change anything below if code structure is not modified.

FILE (GLOB TARGET_SRC "../*.cc", "${srcDefect_DIR}/*.cc", "${src_DIR}/model/dislocation/*.cc" "${src_DIR}/base/IGA_dislocation.cc" )
##---------------------------------------------------------------------------------------------------------------------


DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()


#Include all pathes for head files and link all libraries, modify when you include external libs differently
##---------------------------------------------------------------------------------------------------------------------
INCLUDE_DIRECTORIES(${IGA_INC} ${src_INC} ${hdfInclude} ${SUPERLU_INC} ${szip_INC})
TARGET_LINK_LIBRARIES(${TARGET} ${IGA_LIB} ${hdfLibrary} ${hdfLibrary2} ${SUPERLU_LIB} ${szip_LIB})
