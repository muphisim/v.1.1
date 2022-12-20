This program has been developed on UNIX system. 

Building the source code requires mpi, make, gcc, g++, gfortran, lapack, blas, metis, and petsc. When considering FSI, preCICE is required. 

--------------------------------------------
*Build the code from command line
--------------------------------------------

In the root folder

- without FSI: type "make -j np TYPE=SEQUENTIAL/PARALLEL" in the root directory to obtain the executable file in the bin folder.
- with FSI: type "make -j np TYPE=SEQUENTIAL/PARALLEL PRECICE_DIR=path/to/precice-include" in the root directory to obtain the executable file in the bin folder.
- The code can be compiled without METIS using "METIS=NO" is the runtime parallelisation is not used.
- The code can be compiled with meshless method using "MM=YES". The Sukumar's maxent library needs to be provided by default in path-to-MuPhiSim/lib or by using "MAXENT_LIB_PATH=path/to/maxent-lib" to specify its path.

when TYPE=SEQUENTIAL is used, PETSCMPI=YES must be used to compile the code if pestc is complied in parallel

Once the code is compile sucessfully, update PATH to make the executable MuPhiSim available anywhere by adding the following line to .bashrc

export PATH=path-to-MuPhiSim/bin:$PATH

--------------------------------------------
*Testing
--------------------------------------------

- in sequential: make test TYPE=SEQUENTIAL
- in parallel: make test TYPE=PARALLEL NCPUS=5

--------------------------------------------
*Build documentation
--------------------------------------------

make documentation

View path-to-MuPhiSim/doc/html/index.html to see the documentation.

--------------------------------------------
*Clean
--------------------------------------------

make clean

All existing compilation will be cleaned.

--------------------------------------------
*Run a simulation
--------------------------------------------

Note that the two folders (named by input and output) must be created in the current folder and the input file must be included in  the "input" folder.

- in sequential

MuPhiSim -input example1.inp

- in parallel

mpirun -n 4 MuPhiSim -input example2.inp


