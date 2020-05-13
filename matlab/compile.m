% This is for MATLAB Compiler Runtime (MCR).
mkdir exec
if verLessThan('matlab', '8.5')
    mcc -m OuterHull.m -d exec -a VOXELISE -a ICP_finite -R -nodisplay -R -nojvm
else
    mcc -m OuterHull.m -d exec -a VOXELISE -R -nodisplay -R -nojvm
end
movefile exec/OuterHull ../script
delete exec/*
rmdir exec
