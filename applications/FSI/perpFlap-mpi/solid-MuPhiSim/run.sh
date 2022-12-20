#!/bin/sh

rm -rf output && mkdir output
mpiexec -np 4 ../../../../bin/MuPhiSim-mpi -input perpFlap.inp
