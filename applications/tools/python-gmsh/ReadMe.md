##To add the python path and path to .bashrc

```
./addPythonPath.sh
```


##To view the mesh in a *.inp file using gmsh

```
togmsh.py file.inp
```

##To partion the current inp with N processors using Gmsh
```
partition.py file.inp file-new.inp N
```
where ```file-new.inp``` is the name of new input files with N groups of nodes and elements for parallel simulation with N processors.

## pre-process using Gmsh and Python
- Prepare Model
- Prepare configure

see examples for more information

