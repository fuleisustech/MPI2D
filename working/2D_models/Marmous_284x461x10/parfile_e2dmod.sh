#model_dir=NM_101x301x1
#MODEL PARAMETER
NX=461
NZ=284
DX=10
DZ=10

# VP/VS FILE ONLY USED IN ELASTIC MODELING
VP_FILE= /project/k1046/ful/PhaInv_SeisF90/model/2D_models/Marmous_284x461x10/vel.bin
VS_FILE= /project/k1046/ful/PhaInv_SeisF90/model/2D_models/Marmous_284x461x10/vs.bin
DEN_FILE= /project/k1046/ful/PhaInv_SeisF90/model/2D_models/Marmous_284x461x10/den.bin

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
COORD_IN_FILE= /project/k1046/ful/PhaInv_SeisF90/model/2D_models/Marmous_284x461x10/coord_csg.dat
#DATA_TYPE=w
DATA_TYPE=p

#SOURCE PARAMETER
IsSource=.false.
SOURCE_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/a2d_csg/source.bin
FREQ=15
NW=4000
DTW=0.001
SOURCE_TYPE=P

#FD PARAMETER
FD_TYPE=2ND
FD_ORDER=24
BC_TYPE=1
NPML=400
IsFS=.false.
IsDipole=.false.
IsWriteData=.false.

#RESULT PARAMETER
LOG_FILE=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/e2d_csg/e2dmod.log
#CSG_OUT_W=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/e2d_csg/csg_w_
#CSG_OUT_U=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/e2d_csg/csg_u_
CSG_OUT_P=/project/k1046/ful/PhaInv_SeisF90/results/2D_models/Marmous_284x461x10/e2d_csg/csg_p_

#PARALLEL PARAMETER
IsShotConti=.true.


