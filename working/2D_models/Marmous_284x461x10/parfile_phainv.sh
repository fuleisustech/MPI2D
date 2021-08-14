####################################################
# ITERATION PARAMETER
NIT=200       # maximum iterations

NMAX=2  # max integrations num
TOLER= 0.01   # adjacent misfit change 

NF=6     # num. of high frequencies

FL1=0.1 # 2nd high frequency
FL2=0.1 # 3rd high frequency
FL3=0.1 # 4th high frequency
FL4=0.1 # 4th high frequency
FL5=0.1 # 4th high frequency
FL6=0.1 # 4th high frequency

FH1=3.0 # 2nd high frequency
FH2=5.0 # 3rd high frequency
FH3=7.0 # 4th high frequency
FH4=9.0 # 4th high frequency
FH5=14.0 # 4th high frequency
FH6=20.0 # 4th high frequency

#MODEL PARAMETER
NZ=284
NX=461
DX=10

VEL_FILE=/home/ful/project/MPI2D/results/2D_models/Marmous_284x461x10/vp_ini.bin

IsMask=.false.
MASK_FILE=/home/ful/project/MPI2D/model/2D_models/Marmous_284x461x10/bottom_mask.bin

VMIN=1500
VMAX=4700

IsSeisRecons=.false.
CSG_RECONS=/home/ful/project/MPI2D/results/2D_models/Marmous_284x461x10/mpi/csg_recon_

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
COORD_IN_FILE=/home/ful/project/MPI2D/model/2D_models/Marmous_284x461x10/coord_csg_full.dat
CSG_OUT=/home/ful/project/MPI2D/results/2D_models/Marmous_284x461x10/a2d_csg/csg_
DATA_TYPE=p

IsTaper_dir=.false.
N0=30
T_P=260

IsCut=.false.
NEAR_CUT=5
FAR_CUT=400

#SOURCE PARAMETER
IsSource=.false.
SOURCE_FILE=/home/ful/project/MPI2D/results/2D_models/Marmous_284x461x10/a2d_csg/source.bin
FREQ=15
NW=4000
DTW=0.001

#IS_CONV_SOU=1
#GRAD_FLIP=.true.

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
LOG_FILE=/home/ful/project/MPI2D/results/2D_models/Marmous_284x461x10/mpi/phainv.log
IsPreStackImg=.false.
PRE_IMG_OUT=/home/ful/project/MPI2D/results/2D_models/Marmous_284x461x10/mpi/img_is_
IMAGE_FILE=/home/ful/project/MPI2D/results/2D_models/Marmous_284x461x10/mpi/vel_
GREEN_OUT=/home/ful/project/MPI2D/results/2D_models/Marmous_284x461x10/mpi/g_
GK1_FILE=/home/ful/project/MPI2D/results/2D_models/Marmous_284x461x10/mpi/gk_
# RESIDUAL FILE
RES_FILE=/home/ful/project/MPI2D/results/2D_models/Marmous_284x461x10/mpi/phainv_res.bin

IMAGE_ILLUM_FILE=/home/ful/project/MPI2D/results/2D_models/Marmous_284x461x10/mpi/illum_

#PARALLEL PARAMETER
IsShotConti=.true.
