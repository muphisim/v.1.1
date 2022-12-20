#!/bin/bash

curDir=$PWD
echo $curDir
python model-quad.py

Run=1
if [ $Run -eq 1 ]; then
rm -rf $curDir/muphisim-quad  
mkdir $curDir/muphisim-quad
cd $curDir/muphisim-quad && mkdir input output
cd $curDir && cp input.inp $curDir/muphisim-quad/input
cd $curDir/muphisim-quad && MuPhiSim -input input.inp
echo "done simulation in muphisim"
fi

if [ $Run -eq 1 ]; then
rm -rf $curDir/muphisim-mpi-quad
mkdir $curDir/muphisim-mpi-quad
cd $curDir/muphisim-mpi-quad && mkdir input output 
cd $curDir && cp input.inp $curDir/muphisim-mpi-quad/input
cd $curDir/muphisim-mpi-quad && mpiexec -np 4 MuPhiSim-mpi -input input.inp
echo "done simulation in muphisim-mpi"
fi
