#! /usr/bin/env bash
if [ -d data ]; then
    echo hello
    rm data/*.dat
    rm data/*.csv
fi

dt=1e-15
steps=1e6
J=2
lambda=1e-3
B=10
anisotropyAxis=0
anisotropyPlane=0
T=0
init=2
angle=0
mode=2
structure=chain
nCellsX=30
periodicBoundary=true
stabilize=false

./model.out -dt $dt -steps $steps -J $J -lambda $lambda -B $B \
-anisotropyAxis $anisotropyAxis -anisotropyPlane $anisotropyPlane \
-T $T -init $init -angle $angle -mode $mode -nCellsX $nCellsX \
-structure $structure -periodicBoundary $periodicBoundary -stabilize $stabilize
