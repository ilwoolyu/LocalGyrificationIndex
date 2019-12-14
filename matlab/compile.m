% This is for MATLAB Compiler Runtime (MCR).
% MCR packages are required: MCR_R2015b_glnxa64_installer.zip
mkdir exec
mcc -m OuterHull.m -d exec -a VOXELISE
movefile exec/OuterHull ../script
delete exec/*
rmdir exec
