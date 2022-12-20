
- **examples**: different numerical simulations    
    - create a folder with an arbitrary name
    - change directory to this folder, create two folders: **input** and **output**
    - copy an input file from **examples** to **input**  
    - in terminal, run  
        - in sequential:  
        
            ```muphisim-executable -input input-file.inp```
        
        - in parallel with 8 processors: 
        
            ```mpiexec -np 8 muphisim-executable -input input-file.inp```
    
    - all results are stored in the current folder and **output**
        

- **FSI**: Fluid-Struction-Interaction examples in which the muphisim is used for solid simulation and OpenFOAM is used for fluid simulation. The coupling is performed in the preCICE environment.

- **AM**: 3D printing example.

- **tools**
    - **GUI-Abaqus**: preprocessing using Abaqus through GUI (FSI, stochastic, extradof, etc. not supported)
    - **Hex2Tet-Parallel**: hexahedra-to-tetrahedra conversion
    - **python-gmsh**: pre-processing using Gmsh and python
