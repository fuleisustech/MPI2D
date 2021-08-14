####################################################
# ITERATION PARAMETER
NIT=100       # maximum iterations

NMAX=1  # max integrations num
TOLER= 0.02   # adjacent misfit change 

NF=1     # num. of high frequencies

FL1=0.1 # 2nd high frequency
FL2=4.0 # 3rd high frequency
FL3=6.0 # 4th high frequency
FL4=8.0 # 4th high frequency
FL5=12.0 # 4th high frequency
FL6=12.0 # 4th high frequency

FH1=60.0 # 2nd high frequency
FH2=9.0 # 3rd high frequency
FH3=13.0 # 4th high frequency
FH4=19.0 # 4th high frequency
FH5=50.0 # 4th high frequency
FH6=50.0 # 4th high frequency

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
CSG_RECONS=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/mpi_hil_12_fine/csg_recon_

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
CSG_TMP=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/mpi_hil_12_fine/csgtmp_
DATA_TYPE=p

IsTaper_dir=.true.
N0=30
T_P=260

IsCut=.true.
NEAR_CUT=0
FAR_CUT=200
NO_OFSET=12
OFSET_STEP=400
MAX_OFSET=4600

#SOURCE PARAMETER
IsSource=.false.
SOURCE_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/a2d_csg/source.bin
FREQ=15
NW=4000
DTW=0.001

#FD PARAMETER
FD_TYPE=2ND
FD_ORDER=24
BC_TYPE=1
NPML=200

IsFS=.false.
IsDipole=.false.
IsReadData=.false.

#RTM PARAMETER
IsSmooth_gk=.false.
HW=3

IsDepthPC=.false.
IsIllum=.true.
Illum_Order=2
IsSaveBC=.true.
IsSaveWF=.false.

#RESULT PARAMETER
LOG_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/mpi_hil_12_fine/phainv.log
IsPreStackImg=.false.
PRE_IMG_OUT=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/mpi_hil_12_fine/img_is_
IMAGE_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/mpi_hil_12_fine/vel_
GREEN_OUT=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/mpi_hil_12_fine/g_
GK1_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/mpi_hil_12_fine/gk_
# RESIDUAL FILE
RES_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/mpi_hil_12_fine/phainv_res.bin

IMAGE_ILLUM_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/mpi_hil_12_fine/illum_

#PARALLEL PARAMETER
IsShotConti=.true.
