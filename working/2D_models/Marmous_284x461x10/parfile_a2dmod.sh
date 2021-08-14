####################################################
#MODEL PARAMETER
NX=461
NZ=284
DX=10
DZ=10
VEL_FILE=/home/ful/Project/MPI2D/model/2D_models/Marmous_284x461x10/vel.bin
#VEL_FILE=/home/ful/Project/MPI2D/model/2D_models/Marmous_284x461x10/vp_homo.bin
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
COORD_IN_FILE=/home/ful/Project/MPI2D/model/2D_models/Marmous_284x461x10/coord_csg_full.dat

#SOURCE PARAMETER
#IsSource=.false.
IsSource=.true.
SOURCE_FILE=/home/ful/Project/MPI2D/results/2D_models/Marmous_284x461x10/a2d_csg/source.bin
FREQ=15
NW=4000
DTW=0.001

#FD PARAMETER
FD_TYPE=2ND
FD_ORDER=28
BC_TYPE=1
NPML=200

IsFS=.false.
IsDipole=.false.
IsWriteData=.false.

#RESULT PARAMETER
LOG_FILE=/home/ful/Project/MPI2D/results/2D_models/Marmous_284x461x10/a2d_csg/a2dmod.log
CSG_OUT=/home/ful/Project/MPI2D/results/2D_models/Marmous_284x461x10/a2d_csg/csg_

#PARALLEL PARAMETER
IsShotConti=.true.
