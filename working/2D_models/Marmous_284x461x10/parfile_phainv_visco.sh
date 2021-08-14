####################################################
# ITERATION PARAMETER
NIT=250        # maximum iterations

NMAX=3  # max integrations num
TOLER= 0.015   # adjacent misfit change 

NF=6     # num. of high frequencies

FL1=0.1  # 1st high frequency
FL2=2.0 # 2nd high frequency
FL3=5.00 # 3rd high frequency
FL4=9.0 # 4th high frequency
FL5=14.0 # 4th high frequency
FL6=20.0 # 4th high frequency

FH1=5.0  # 1st high frequency
FH2=9.0 # 2nd high frequency
FH3=15.0 # 3rd high frequency
FH4=19.0 # 4th high frequency
FH5=24.0 # 4th high frequency
FH6=50.0 # 4th high frequency

#MODEL PARAMETER
NZ=284
NX=461
DX=10

VEL_FILE=/project/k1046/ful/PhaInv_SeisF90/model/2D_models/Marmous_284x461x10/vsmooth_fine.bin


IsMask=.true.
#MASK_FILE=/project/k1046/ful/PhaInv_SeisF90/model/2D_models/Marmous_284x461x10/mask_bottom.bin
MASK_FILE=/project/k1046/ful/PhaInv_SeisF90/model/2D_models/Marmous_284x461x10/mask_bottom_e.bin

VMIN=1500
VMAX=4700

#IsSeisRecons=.true.
IsSeisRecons=.false.
CSG_RECONS=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/phainv_visco2/csg_recon_

#INV_TYPE=CG
INV_TYPE=SD

#DEN_FILE=/home/guob/projects/SeisF90/model/2D_models/salt_150fx645fx30f/den.bin
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
COORD_IN_FILE=/project/k1046/ful/PhaInv_SeisF90/model/2D_models/Marmous_284x461x10/coord_csg.dat

CSG_OUT=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/a2d_csg_visco2/csg_

CSG_TMP=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/phainv_visco2/csgtmp_
DATA_TYPE=p

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
NPML=400
IsFS=.false.
IsDipole=.false.
IsReadData=.false.

#RTM PARAMETER
IsIllum=.true.
Illum_Order=2
IsSaveBC=.true.
IsSaveWF=.false.


IsSmooth_gk=.true.
HW=1

#RESULT PARAMETER
LOG_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/phainv_visco2/phainv.log
IsPreStackImg=.false.
PRE_IMG_OUT=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/phainv_visco2/img_is_
IMAGE_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/phainv_visco2/vel_
GREEN_OUT=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/phainv_visco2/g_
GK1_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/phainv_visco2/gk_

# RESIDUAL FILE
RES_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/phainv_visco2/phainv_res.bin
IMAGE_ILLUM_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/phainv_visco2/illum_

#PARALLEL PARAMETER
IsShotConti=.true.

