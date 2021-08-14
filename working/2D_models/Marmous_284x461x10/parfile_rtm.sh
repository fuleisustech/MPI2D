####################################################
#MODEL PARAMETER
NZ=284
NX=461
DX=10

NH=41

VEL_FILE=/project/k1046/ful/PhaInv_SeisF90/model/2D_models/Marmous_284x461x10/vel.bin
#VEL_FILE=/project/k1046/ful/PhaInv_SeisF90/model/2D_models/Marmous_284x461x10/vsmooth_bad.bin
#VEL_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/fwi_ms_bad/vel_24.bin
#VEL_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/phainv_bad/vel_145.bin

IsAper=.false.
X0=0
Z0=0
NT=4000
NT_OUT=4000
DT=0.001
DT_OUT=0.001
SORT_IN=CSG
COORDFILE_TYPE=ASCII
COORD_IN_FILE=/project/k1046/ful/PhaInv_SeisF90/model/2D_models/Marmous_284x461x10/coord_csg.dat

CSG_OUT=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/a2d_csg/csg_

#SOURCE PARAMETER
IsSource=.false.
SOURCE_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/a2d_csg/source.bin
FREQ=20
NW=4000
DTW=0.001

#FD PARAMETER
FD_TYPE=2ND
FD_ORDER=28
BC_TYPE=1
NPML=400
IsFS=.false.
IsDipole=.false.
IsReadData=.false.
#IMAGE_CONDITION=2
IMAGE_CONDITION=0

#RTM PARAMETER
IsIllum=.true.
#IsIllum=.false.
IsSaveBC=.true.
IsSaveWF=.false.

#RESULT PARAMETER

LOG_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/rtm_h/rtm.log

#IsPreStackImg=.true.
IsPreStackImg=.false.
PRE_ILLUM_IMG_OUT=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/rtm_h/img_is_
PRE_IMG_OUT=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/rtm_h/img_is_


IMAGE_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/rtm_h/img_h_true.bin
#IMAGE_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/rtm_h/img_h_ini.bin
#IMAGE_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/rtm_h/img_h_fwi.bin
#IMAGE_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/rtm_h/img_h_mpi.bin

#PARALLEL PARAMETER
IsShotConti=.true.

