####################################################
# ITERATION PARAMETER
NIT=80

TOLER= 0.01

#MODEL PARAMETER
NZ=284
NX=461
DX=10

VEL_FILE=/project/k1046/ful/PhaInv_SeisF90/model/2D_models/Marmous_284x461x10/vp0_fine.bin
#VEL_FILE=/project/k1046/ful/PhaInv_SeisF90/model/2D_models/Marmous_284x461x10/vp0_bad.bin

IsMask=.true.
MASK_FILE=/project/k1046/ful/PhaInv_SeisF90/model/2D_models/Marmous_284x461x10/bottom_mask.bin

VMIN=1500
VMAX=4700

IsSeisRecons=.false.
CSG_RECONS=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/fwi_hil_zedong12_fine_fs/csg_recon_

IsRemove_dir=.false.
#INV_TYPE=CG
INV_TYPE=SD

#Data PARAMETER
IsAper=.false.
X0=0
Z0=0
NT=4000
NT_OUT=4000
DT=0.001
DT_OUT=0.001

SORT_IN=CSG
COORDFILE_TYPE=ASCII
COORD_IN_FILE=/project/k1046/ful/PhaInv_SeisF90/model/2D_models/Marmous_284x461x10/coord_csg_full.dat
CSG_OUT=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/a2d_csg/csg_
CSG_DIR=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/a2d_csg/csg0_
CSG_TMP=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/fwi_hil_zedong12_fine_fs/csgtmp_
DATA_TYPE=p

IsTaper_dir=.true.
N0=30
T_P=300

IsCut=.true.
NEAR_CUT=5
FAR_CUT=400

#SOURCE PARAMETER
IsSource=.false.
SOURCE_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/a2d_csg_fs/source.bin
FREQ=15
NW=4000
DTW=0.001

#FD PARAMETER
FD_TYPE=2ND
FD_ORDER=24
BC_TYPE=1
NPML=400

IsFS=.true.
IsDipole=.false.
IsReadData=.false.

#RTM PARAMETER
IsIllum=.true.
Illum_Order=2
IsSaveBC=.true.
IsSaveWF=.false.

IsSmooth_gk=.false.
HW=3

#RESULT PARAMETER
LOG_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/fwi_hil_zedong12_fine_fs/fwi.log
IsPreStackImg=.false.
PRE_IMG_OUT=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/fwi_hil_zedong12_fine_fs/img_is_
IMAGE_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/fwi_hil_zedong12_fine_fs/vel_
GREEN_OUT=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/fwi_hil_zedong12_fine_fs/g_
GK1_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/fwi_hil_zedong12_fine_fs/gk_

# RESIDUAL FILE
RES_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/fwi_hil_zedong12_fine_fs/fwi_res.bin
IMAGE_ILLUM_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/fwi_hil_zedong12_fine_fs/illum_

#PARALLEL PARAMETER
IsShotConti=.true.


