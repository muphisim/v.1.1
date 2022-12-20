## File CreateInputFile.py

The main file to create input file for MuPhiSim. It is mandatory to provide the config (a Configure instance) and model (a Model instance).

Run **addPythonPath.sh** in terminal to make current folder available in $PATH and $PYTHONPATH

## File Configure.py

This is a virtual class to define a configuration including the FEM configuration, the solver, the boundary conditions, material laws, etc. 

## File Model.py

This is a virtual class to define a finite element mesh. See comments in Model.py for more details.

A derived class based on Gmsh **GmshModel.py** was build.
 
## File partion.py

To partition the inp file using 

**```partion.py input.inp output.inp N```**

where N is the number of partitions.  

## File togmsh.py

To view the inp file using ```togmsh your-inp-file.inp``` in the terminal.
