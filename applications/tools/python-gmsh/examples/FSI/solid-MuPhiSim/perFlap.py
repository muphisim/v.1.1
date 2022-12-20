import gmsh
import GmshModel
import Configure
import CreateInputFile


def PerFlap():
    # geometrical parameters
    H = 1
    W = 0.1
    n_x_Direction = 10
    n_y_Direction = 30
            
    model = gmsh.model
    model.add("PerFlap")
    
    lc = 1e-1*W
    model.geo.addPoint(-0.5*W, 0, 0, lc, 1)
    model.geo.addPoint(0.5*W, 0, 0, lc, 2)
    model.geo.addPoint(0.5*W, H, 0, lc, 3)
    model.geo.addPoint(-0.5*W, H, 0, lc, 4)
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
    
    model.geo.addPhysicalGroup(1, [1], 1)
    model.geo.addPhysicalGroup(1, [2,3,4], 2)
    model.geo.addPhysicalGroup(2, [1], 10)
    
    model.geo.synchronize()

    model.mesh.generate(model.getDimension())
    model.mesh.setOrder(2) 
    return model

class myConfigure (Configure.Configure):
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
                #scale factor to estimate time sep for solver
                "Scale Factor": 1e10,
                "FSI": ((1,2),),
                #output specification, each component of the list corresponds a line in the oxfemm input file
                "Outputs": ((1000,), ("All" , "max"),),
                # start and end time
                "Time": (0.,5.),
                "DispBC":(
                            # each bc contains (type, location(dim, physical) list of prescried components and prescribed values
                            ("DISPLACEMENT RAMP",(1,1), ((0,0.),(1,0.))), 
                         ), 
                },
                )
        # define all Neumann BCs
        self.NeumannBCs = {
                "PressureBC": (
                                # type, location, time interval, value
                                #("PRESSURE RAMP", (1,2), (0.,10.), -1.),
                                #("PRESSURE RAMP", (2,52), (0.5,1.), -0.1),
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
                "Density": 3000.,
                #material type
                "Type": "Linear-Elastic",
                # parameter of the law, each component of this list corresponds a line in the input file
                "Parameters": [(4e6,0.3,"PlaneStrain"),],
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
model = PerFlap()
gmsh.write("perFlap.msh")
gmshModel = GmshModel.GmshModel(model)
config = myConfigure()
CreateInputFile.createInputFile(fileName="perFlap.inp",model=gmshModel,configure=config)
gmsh.finalize()



