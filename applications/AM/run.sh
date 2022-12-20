rm -rf input output
mkdir input output
cp input.inp input
mpiexec -np 5 MuPhiSim-mpi -input input.inp
