#! /usr/bin/env bash
if [ -d data ]; then
    echo hello
    rm data/*.dat
    rm data/*.csv
fi

dt=1e-15
steps=1e6
J=0
lambda=1e-3
B=0
anisotropyAxis=0
anisotropyPlane=0
T=1
init=1
angle=45
mode=1
structure=single
nCellsX=1
periodicBoundary=true
dipole=true
stabilize=false

./model.out -dt $dt -steps $steps -J $J -lambda $lambda -B $B \
-anisotropyAxis $anisotropyAxis -anisotropyPlane $anisotropyPlane \
-T $T -init $init -angle $angle -mode $mode -nCellsX $nCellsX \
-structure $structure -periodicBoundary $periodicBoundary -dipole $dipole -stabilize $stabilize
#Added .out to correspond to file