import gmsh
import GmshModel
import Configure
import CreateInputFile

def createGmshModel(name):
    # geometrical parameters
    Lx = 1
    Ly = 0.1
    Lz = 0.04
    n_x_Direction = 15
    n_y_Direction = 4
    n_z_Direction = 3
            
    model = gmsh.model
    model.add(name)

    lc = 1e-1*Lx
    model.geo.addPoint(0, 0, 0, lc, 1)
    model.geo.addPoint(Lx, 0, 0, lc, 2)
    model.geo.addPoint(Lx, Ly, 0, lc, 3)
    model.geo.addPoint(0, Ly, 0, lc, 4)
    model.geo.addLine(1, 2, 1)
    model.geo.addLine(3, 2, 2)
    model.geo.addLine(3, 4, 3)
    model.geo.addLine(4, 1, 4)
    model.geo.addCurveLoop([4, 1, -2, 3], 1)
    model.geo.addPlaneSurface([1], 1)

    model.geo.mesh.setTransfiniteCurve(1, n_x_Direction+1)
    model.geo.mesh.setTransfiniteCurve(3, n_x_Direction+1)
    model.geo.mesh.setTransfiniteCurve(2, n_y_Direction+1)
    model.geo.mesh.setTransfiniteCurve(4, n_y_Direction+1)
    model.geo.mesh.setTransfiniteSurface(1, "Right", [1, 2, 3, 4])

    model.geo.addPhysicalGroup(2, [1], 1)
    model.geo.addPhysicalGroup(1, [1], 11)
    model.geo.addPhysicalGroup(1, [2], 12)
    model.geo.addPhysicalGroup(1, [3], 13)
    model.geo.addPhysicalGroup(1, [4], 14)

    model.geo.synchronize()

    model.mesh.generate(model.getDimension())
    #model.mesh.setOrder(2) 
    return model

class myConfigure(Configure.Configure):
    def __init__(self):
        """
        configure in oxfem
        Attributes:
            FEM     : all FEM locations
            MM      : all MM locations
            solvers : all solvers
            materials: all material behaviors
                
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
                "Solver": "IMPLICIT",
                # number of step
                "Number Of Steps": 70,
                #scale factor to estimate time sep for solver
                "Scale Factor": 500,
                #output specification, each component of the list corresponds a line in the oxfemm input file
                "Outputs": ((1000,), ("All" , "max"),),
                # start and end time
                "Time": (0.,4.),
                "DispBC":(
                            # each bc contains (type, location(dim, physical) list of prescried components and prescribed values
                            ("DISPLACEMENT RAMP",(1,14), ((0,0.),(1,0.))), 
                         ), 
                "Extraction":(
                               ("NODE","TipDisplacement","Unknown","Mean",78,0),
                               ("NODE","TipDisplacement","Unknown","Mean",78,1),
                               ("NODE","TipInternalForce","InternalForce","Sum",(1,12),0),
                               ("NODE","TipInternalForce","InternalForce","Sum",(1,12),1),
                               ("NODE","TipExternalForce","ExternalForce","Sum",(1,12),0),
                               ("NODE","TipExternalForce","ExternalForce","Sum",(1,12),1),
                               ("NODE","RootInternalForce","InternalForce","Sum",(1,14),0),
                               ("NODE","RootInternalForce","InternalForce","Sum",(1,14),1),
                               ("NODE","RootExternalForce","ExternalForce","Sum",(1,14),0),
                               ("NODE","RootExternalForce","ExternalForce","Sum",(1,14),1),
                               ("ELEMENT","DeforEnergy","DEFO_ENERGY","Sum","All"),
                               ("ELEMENT","kineEnergy","KIN_ENERGY","Sum","All"),
                               ("ELEMENT","extEnergy","EXTERNAL_ENERGY","Sum","All"),
                            )
                },
                )
        # define all Neumann BCs
        self.NeumannBCs = {
                "PressureBC": (
                                # type, location, time interval, value
                                #("PRESSURE RAMP", (2,13), (0.,1.), 1.),
                              ),
                "TractionBC": (
                                 # type, location, time interval, values
                                ("TRACTION RAMP", (1,12), (0.,0.8001), (0,1.,0)),
                              ),
                "Flux ExtraDof": (
                                    # type, location, time interval, value
                                    #("HEAT INST", (2,52), (0., 0.5), (800.,800.,0.)),
                                    # type, location, time interval, value
                                    #("VOLUMETRIC HEAT FLUX", (3,1), (0.,1.), 100.),
                                    # type, location, time interval, h and sink temperature
                                    #("CONVECTION", (2,52), (0., 0.5), (800.,0)),
                                    # type, location, time interval, A and sink temperature
                                    #("RADIATION", (2,52), (0., 0.5), (500.,0)),
                                  )
                #"GRAVITY": 1.
                }
        # define material laws  - a tuple of dictionary
        self.materials = ({
                #material name
                "Name": "Mat1",
                #density
                "Density": 1.,
                #material type
                "Type": "Linear-Elastic",
                # parameter of the law, each component of this list corresponds a line in the input file
                "Parameters": [(1000,0.3,"PlaneStrain"),],
                # support All or [(dim,phys), ...]
                #"ExtraDof": ("Temperature",(1.e-10, 0, 5, 0, 0, 0, 0)),
                "List of Elements": "All"
                #"List of Elements": [(3,1),]
                },
                )
        self.initialBCs = {
                #"INITIAL CONDITIONS EXTRADOF": (273.,)
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

model = createGmshModel("myModel")
gmsh.write("mymodel.msh")
gmshModel = GmshModel.GmshModel(model)
config = myConfigure()

CreateInputFile.createInputFile(fileName="input.inp",model=gmshModel,configure=config)

gmsh.finalize()



