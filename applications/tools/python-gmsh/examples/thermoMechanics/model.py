import gmsh
import GmshModel
import Configure
import CreateInputFile
import numpy as np

class myConfigure(Configure.Configure):
    def __init__(self):
        """
        configure in oxfem
        Attributes:
            FEM     : all FEM locations
            MM      : all MM locations
            solvers : all solvers
            materials: all material behaviors
            NeumannBCs: all Neumann BCs
            initialBCs: all initial BCs
        """
        # define FEM domains
        self.FEM = {
                   #support All or [(dim,phys), ...] or nothing
                   "Support": "All"
                   #"Support": [(3,1)]
                   }
        # define MM domains
        self.MM = {
              #support All or [(dim,phys), ...] or nothing
              #"Support": "All"
              #"Support": [(3,1)]
              }
        # define solvers - a tuple of dictionary
        self.solvers = ({
                #solver type
                "Solver": "IMPLICIT STATIC",
                #scale factor to estimate time sep for solver
                #"Scale Factor": 1e5,
                "Number of steps":100,
                #output specification, each component of the list corresponds a line in the oxfemm input file
                "Outputs": ((1000,), ("All" , "max"),),
                # start and end time
                "Time": (0.,1e4),
                "DispBC":(
                            # each bc contains (type, location(dim, physical) list of prescried components and prescribed values
                            ("DISPLACEMENT INST",(2,6), ((0,0.),(1,0.))), 
                            ("DISPLACEMENT INST",(1,11), ((2,10.),)), 

                         ), 
                "Extraction":(
                          
                            )
                },
                )
        # define all Neumann BCs
        self.NeumannBCs = {
                "PressureBC": (
                                
                              ),
                "Flux ExtraDof": (
                                    # type, location, time interval, value
                                    #("HEAT INST", (2,52), (0., 5), (0.,0.,0.)),
                                    # type, location, time interval, value
                                    #("VOLUMETRIC HEAT FLUX", (3,80), (0.,10), 1e3),
                                    # type, location, time interval, h and sink temperature
                                    #("CONVECTION", (2,52), (0., 7.5), (800.,0)),
                                    # type, location, time interval, A and sink temperature
                                    #("RADIATION", (2,52), (0., 8.5), (500.,0)),
                                  )
                #"GRAVITY": 1.
                }
        # define material laws  - a tuple of dictionary
        self.materials = ({
                #material name
                "Name": "Mat1",
                #densityG
                "Density": 2700.,
                #material type
                "Type": "Linear-Elastic",
                #"List of Elements": [(3,1),]
                # parameter of the law, each component of this list corresponds a line in the input file
                "Parameters": [(70e3,0.3, "PlaneStrain"),],
                # extra dof field
                "ExtraDof": ("Temperature",(910e-6, 0, 237e-6, 0, 293., 293., 293.)),
                # support All or [(dim,phys), ...]
                "List of Elements": "All",
                },
                )
        self.initialBCs = {
                "INITIAL CONDITIONS EXTRADOF": (293.,)
                }
                
    def getFEMConfigure(self):
        return self.FEM
        
    def getMMConfigure(self):
        return self.MM
    
    def getSolversConfigure(self):
        return self.solvers
        
    def getNeumannBCsConfigure(self):
        return self.NeumannBCs
    
    def getMaterialsConfigure(self):
        return self.materials 
    
    def getInitialBCsConfigure(self):
        return self.initialBCs


gmsh.initialize()
gmsh.open("model.msh")
gmshModel = GmshModel.GmshModel(gmsh.model())
print(gmshModel)
print(gmshModel.getGroupOfNodes((1,11)))
config = myConfigure()
CreateInputFile.createInputFile(fileName="input.inp",model=gmshModel,configure=config)
gmsh.finalize()


