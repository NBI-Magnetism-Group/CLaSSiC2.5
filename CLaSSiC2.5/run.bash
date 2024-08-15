#! /usr/bin/env bash

if [ -d data ]; then
    echo hello
    rm data/*.dat
    rm data/*.csv
fi

dt=1e-17
steps=1e6
J=-0.1 #-2.4
Jz=-1 #-2.4
Jx=-1 #-2.4
lambda=1e-3
B=0
anisotropyAxis=0 #-0.13 #Nicklas: set to meV 
anisotropyPlane=0 #just use negative axis instead.
T=30
init=9
angle=15
mode=2
structure=square
nCellsX=20
periodicBoundary=true
#dipole=false
stabilize=false

./model.out -dt $dt -steps $steps -J $J -Jz $Jz -Jx $Jx -lambda $lambda -B $B \
-anisotropyAxis $anisotropyAxis -anisotropyPlane $anisotropyPlane \
-T $T -init $init -angle $angle -mode $mode -nCellsX $nCellsX \
-structure $structure -periodicBoundary $periodicBoundary -stabilize $stabilize
#Added .out to correspond to file    -dipole $dipole