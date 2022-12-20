import gmsh
import GmshModel
import Configure
import CreateInputFile
import numpy as np

af = 0.055e-3
bf = 0.11e-3
cf = 0.055e-3


P = 180.
eta = 0.6
qmax = 6*np.sqrt(3)*P*eta/(af*bf*cf*np.pi*np.sqrt(np.pi)) # 836113205.51
print(qmax)
x0, y0, z0 = 0.1e-3,0.3e-3,0.2e-3
vx, vy, vz = 0,960e-3,0
a,b,c = af/np.sqrt(3), bf/np.sqrt(3), cf/np.sqrt(3)
print(a,b,c)

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
                "Number of steps":200,
                #output specification, each component of the list corresponds a line in the oxfemm input file
                "Outputs": ((100000,), ("All" , "max"),),
                # start and end time
                "Time": (0.,0.005),
                "DispBC":(
                            # each bc contains (type, location(dim, physical) list of prescried components and prescribed values
                            ("DISPLACEMENT RAMP",(3,80), ((0,0.),(1,0.),(2,0.))), 
                            ("DISPLACEMENT RAMP",(3,81), ((0,0.),(1,0.),(2,0.))), 
                            ("DISPLACEMENT RAMP",(3,82), ((0,0.),(1,0.),(2,0.))), 
                            ("DISPLACEMENT RAMP",(2,86), ((3,0.),)), 
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
                                    #("VOLUMETRIC HEAT FLUX", (3,81), (0.,5e-5), 78288022.),
                                    ("GAUSSIAN VOLUMETRIC HEAT FLUX", "ALL", (0.,1.0417e-03), (qmax,x0, y0, z0,vx, vy,vz,a,b,c)),
                                    # type, location, time interval, h and sink temperature
                                    #("CONVECTION", (2,85), (0., 1), (25.,298)),
                                    # type, location, time interval, A and sink temperature
                                    #("RADIATION", (2,85), (0., 1), (0.5*5.67e-8,298)),
                                  )
                #"GRAVITY": 1.
                }
        # define material laws  - a tuple of dictionary
        self.materials = ({
                #material name
                "Name": "Mat1",
                #densityG
                "Density": 8146.,
                #material type
                "Type": "Linear-Elastic",
                #"List of Elements": [(3,1),]
                # parameter of the law, each component of this list corresponds a line in the input file
                "Parameters": [(200e9,0.3),],
                # extra dof field
                "ExtraDof": ("Temperature",(435, 0, 12, 0, 300., 300., 300.)),
                # support All or [(dim,phys), ...]
                "List of Elements": "All",
                },
                )
        self.initialBCs = {
                "INITIAL CONDITIONS EXTRADOF": (300.,)
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
config = myConfigure()
CreateInputFile.createInputFile(fileName="input.inp",model=gmshModel,configure=config)
gmsh.finalize()


