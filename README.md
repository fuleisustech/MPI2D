# MPI2D
This is a seismic data inversion package for multiscale phase inversion. You may need intel compiler.

# step 1
cd to src/, modify Line 4 according to the location of MPI2D in Makefile.config
MYDIR=/home/ful/Project  -> MYDIR=/the/location/of/MPI2D

cd to src/core/, make clean; make;
cd to src/apps/, make;  

# step 2
cd to model/2D_models/Marmous_284x461x10, 
build the true and initial velocity model, vp_true and vp_ini model; 
and create a coordinate file with mkcoord.m

# step 3
cd working/
make a2dmod, start do forward modeling

