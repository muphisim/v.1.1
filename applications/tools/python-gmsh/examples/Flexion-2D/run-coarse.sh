#!/bin/bash

curDir=$PWD
echo $curDir
python model-coarse.py

Run=1
if [ $Run -eq 1 ]; then
rm -rf $curDir/muphisim-coarse
mkdir $curDir/muphisim-coarse
cd $curDir/muphisim-coarse && mkdir input output
cd $curDir && cp input.inp $curDir/muphisim-coarse/input
cd $curDir/muphisim-coarse && MuPhiSim -input input.inp
echo "done simulation in muphisim"
fi

if [ $Run -eq 1 ]; then
rm -rf $curDir/muphisim-mpi-coarse
mkdir $curDir/muphisim-mpi-coarse
cd $curDir/muphisim-mpi-coarse && mkdir input output 
cd $curDir && cp input.inp $curDir/muphisim-mpi-coarse/input
cd $curDir/muphisim-mpi-coarse && mpiexec -np 4 MuPhiSim-mpi -input input.inp
echo "done simulation in muphisim-mpi"
fi
