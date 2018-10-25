% This is for MATLAB Compiler Runtime (MCR).
mcc -m voxelizer.m -d exec -a VOXELISE
delete exec/*.txt exec/*.sh exec/*.log