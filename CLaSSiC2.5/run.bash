#! /usr/bin/env bash
if [ -d data ]; then
    echo hello
    rm data/*.dat
    rm data/*.csv
fi

dt=1e-15
steps=1e6
J=0
lambda=0
B=5
anisotropyAxis=0
anisotropyPlane=0
T=0
init=1
angle=45
mode=0
structure=single
nCellsX=1
periodicBoundary=true
stabilize=false

./model.out -dt $dt -steps $steps -J $J -lambda $lambda -B $B \
-anisotropyAxis $anisotropyAxis -anisotropyPlane $anisotropyPlane \
-T $T -init $init -angle $angle -mode $mode -nCellsX $nCellsX \
-structure $structure -periodicBoundary $periodicBoundary -stabilize $stabilize
