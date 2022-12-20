
These are external libraries used by MuPhiSim:

- MAXENT-V1.4 folder

The source code was downloaded directly from https://imechanica.org/files/maxentcode.tar_.doc

It contains the Fortran 90 source code for the max entropy shape functions of Sukumar (Sukumar and Wright, 2007).


To compile the code, 
- gfortran is used as compiler by modifying MAXENT-V1.4/src/makefile.inc when specifying G95 as

G95      = gfortran  

- In MAXENT-V1.4/src/makefile, /usr/lib/libg2c.so.0 is not used, so

 EXT_LIBS  = -lblas -llapack 


- In current folder and from terminal, using

mkdir -p lib
make maxentlib

