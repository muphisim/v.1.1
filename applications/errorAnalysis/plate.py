import gmsh
import GmshModel
import Configure
import CreateInputFile
import os
from argparse import ArgumentParser

import numpy as np
import matplotlib.pyplot as plt

E = 2e3
nu = 0.3
c =1.
l = 3.
I = 2.*c*c*c/3.
P = 150.

def getUV(x, y):
    u = (P*y/(6*E*I))*((6*l-3*x)*x+(2+nu)*(y*y-c*c))
    v = -(P/(6*E*I))*(3*nu*y*y*(l-x)+(4+5*nu)*c*c*x+(3*l-x)*x*x)
    return u, v

plt.figure()

x = np.linspace(l,l,1000)
y = np.linspace(-c,c,1000)
u = y*0
v = y*0
for ii in range(1000):
    u[ii], v[ii] = getUV(x[ii], y[ii])
    
plt.plot(x,y,label="initial")
plt.plot(x+u,y+v,label="deformed")

x = np.linspace(0,l,1000)
y = np.linspace(c,c,1000)
u = y*0
v = y*0
for ii in range(1000):
    u[ii], v[ii] = getUV(x[ii], y[ii])
    
plt.plot(x,y,label="initial")
plt.plot(x+u,y+v,label="deformed")

x = np.linspace(0,l,1000)
y = np.linspace(-c,-c,1000)
u = y*0
v = y*0
for ii in range(1000):
    u[ii], v[ii] = getUV(x[ii], y[ii])
    
plt.plot(x,y,label="initial")
plt.plot(x+u,y+v,label="deformed")

x = np.linspace(0,0,1000)
y = np.linspace(-c,c,1000)
u = y*0
v = y*0
for ii in range(1000):
    u[ii], v[ii] = getUV(x[ii], y[ii])
    
plt.plot(x,y,label="initial")
plt.plot(x+u,y+v,label="deformed")
plt.legend()
#plt.show()


class myConfigure(Configure.Configure):
    def __init__(self, bcNodes=[], FEM=[], MM=[]):
        """
        configure in oxfem
        Attributes:
            FEM     : all FEM locations
            MM      : all MM locations
            solvers : all solvers
            materials: all material behaviors
                
        """
        # define FEM domains
        if MM is None:
            self.FEM = {
                       #support All or [(dim,phys), ...] or nothing
                       "Support": "All"
                       #"Support": [(2,ii) for ii in FEM]
                       }
        elif FEM is not None:
            self.FEM = {
                       #support All or [(dim,phys), ...] or nothing
                       #"Support": "All"
                        "Support": [(2,ii) for ii in FEM]
                       }
        else:
            self.FEM ={}
        # define MM domains
        if FEM is None:
            self.MM = {
                  #support All or [(dim,phys), ...] or nothing
                  "Support": "All"
                  #"Support": [(2,ii) for ii in MM]
                  }
        elif MM is not None:
            self.MM = {
                  #support All or [(dim,phys), ...] or nothing
                  #"Support": "All"
                  "Support": [(2,ii) for ii in MM]
                  }
        else:
            self.MM = {}
        # define solvers - a tuple of dictionary
        self.solvers = ({
                #solver type
                "Solver": "IMPLICIT STATIC",
                # number of step
                "Number Of Steps": 1,
                #output specification, each component of the list corresponds a line in the oxfemm input file
                "Outputs": ((1000,), ("All" , 10),),
                # start and end time
                "Time": (0.,1.),
                "Adaptive Time step": ("POWER",8,0.5),
                "Max Iterations": 10,
                "Absolute Tolerance": 1e-8,
                "Relative Tolerance": 1e-6,
                "PETSC_SOLVER":"cg",
                "PETSC_PC":"gamg",
                "ERROR ANALYSIS": ["CantileverBeamUnderParabolicTraction",l,c,E,nu,P],
                "DispBC":(
                            # each bc contains (type, location(dim, physical) list of prescried components and prescribed values
                            #("DISPLACEMENT RAMP",(0,32), ((0,0.),(1,0.))), 
                            #("DISPLACEMENT RAMP",(1,30), ((0,0.),)), 
                         ), 
                 
                "NodeDispBC": [["DISPLACEMENT RAMP",nn[0], ((0,nn[1]),(1,nn[2]))] for nn in bcNodes],
                "Extraction":(
                               ("NODE","forceGr12","InternalForce","Sum",(1,12),1),
                               ("NODE","forceGr12","InternalForce","Sum",(1,12),0),
                               ("NODE","forceGr13","InternalForce","Sum",(1,13),1),
                               ("NODE","forceGr13","InternalForce","Sum",(1,13),0),
                               ("NODE","Group10","Unknown","Rough",(2,10),0),
                               ("NODE","Group10","Unknown","Rough",(2,10),1),
                               ("NODE","Group11","Unknown","Rough",(2,11),0),
                               ("NODE","Group11","Unknown","Rough",(2,11),1),
                            ), 
                },
                )
        # define all Neumann BCs
        self.NeumannBCs = {
                "PressureBC": (
                                # type, location, time interval, value
                                #("PRESSURE RAMP", (2,21), (0.,1.), 1.e1),
                              ),
                "TractionBC": (
                                 # type, location, time interval, values
                                #("TRACTION RAMP", (1,15), (0.,1), (0,-1.e3,0)),
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
                "Density": 1.e-9,
                #material type
                "Type": "Linear-Elastic",
                # parameter of the law, each component of this list corresponds a line in the input file
                "Parameters": [(E,nu,), ],
                # support All or [(dim,phys), ...]
                "List of Elements":"All"
                #"ExtraDof": ("Temperature",(1.e-10, 0, 5, 0, 0, 0, 0)),
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


parser = ArgumentParser()
parser.add_argument("-m", "--mesh", type=str, required=True,
                    help="mesh file name", metavar="MESH")
parser.add_argument("-MMScalingFactor", "--MMScalingFactor", nargs='?',  type=float,
                    help="MM scaling factor", metavar="factor",default=1.)

parser.add_argument("-outFolder", "--outFolder", type=str, required=True,
                    help="output folder", metavar="factor")
                    
parser.add_argument('-FEM','--FEM', nargs="+", type=int)
parser.add_argument('-MM','--MM', nargs="+", type=int)

parser.add_argument("-remesh", "--remesh", nargs='?', type=int,
                    help="remesh option", metavar="remesh",default=0)

parser.add_argument("-remeshBound", "--remeshBound", nargs='?', default=0, type=int,
                    help="remesh option", metavar="remeshBound")

parser.add_argument("-MMIntegOrder", "--MMIntegOrder", nargs='?', default=0, type=int,
                    help="MMIntegOrder", metavar="MMIntegOrder")
                   
parser.add_argument("-MMAdaptiveSupportRadius", "--MMAdaptiveSupportRadius", nargs='?', default=1, type=int,
                    help="MMAdaptiveSupportRadius", metavar="n")

args = parser.parse_args()

print(args)
print(args.remesh)


def getIntersectionElements(model, grLeft, grRight, includeLeft):
    allNode1= model.getGroupOfNodes(grLeft)
    allNode2= model.getGroupOfNodes(grRight)
    commonNodes = set(allNode1)&set(allNode2)
    allElems1 = model.getGroupOfElements(grLeft)
    allElems2 = model.getGroupOfElements(grRight)
    
    allElems=[]
    if includeLeft:
        allElems =allElems+ [i for i in allElems1]
    allElems = allElems+[i for i in allElems2]
    
    refinedElements=[]
    for ele in allElems:
        nodeList = set(model.elements[ele-1])
        if  len(commonNodes&nodeList) >0:
            refinedElements.append(ele)
    
    return refinedElements
   

gmsh.initialize()
gmsh.open(args.mesh)
model =gmsh.model()
gmshModel = GmshModel.GmshModel(model)

if args.remesh == 1:
    allRefineElements= gmshModel.getGroupOfElements((2,11))
    gmshModel.centroid_refine(boundaryOnly=False, elementGroups=allRefineElements)
if args.remeshBound==1:
    #allRefineElements= gmshModel.getGroupOfElements((2,11))
    allRefineElements= getIntersectionElements(gmshModel,(2,10),(2,11),False)
    gmshModel.centroid_refine(boundaryOnly=True, elementGroups=allRefineElements)


bcNodes = []
allNodes = gmshModel.getGroupOfNodes((1,12))
for nn in allNodes:
    coords = gmshModel.nodes[nn-1]
    u, v = getUV(coords[0],coords[1])
    bcNodes.append([nn,u,v])
allNodes = gmshModel.getGroupOfNodes((1,13))
for nn in allNodes:
    coords = gmshModel.nodes[nn-1]
    u, v = getUV(coords[0],coords[1])
    bcNodes.append([nn,u,v])
config = myConfigure(bcNodes=bcNodes, FEM=args.FEM, MM=args.MM)

os.system(f"rm -rf {args.outFolder}")
os.makedirs(args.outFolder, exist_ok = True)

nameSplit = args.mesh.split(".")
inputFileName = nameSplit[0]+".inp"
CreateInputFile.createInputFile(fileName=os.path.join(args.outFolder,inputFileName),model=gmshModel,configure=config)

os.system(f"MuPhiSim -input {inputFileName} -inputDir {args.outFolder} -outputDir {args.outFolder} -MMScalingFactor {args.MMScalingFactor} -MMIntegOrder {args.MMIntegOrder} -MMAdaptiveSupportRadius {args.MMAdaptiveSupportRadius} | tee {args.outFolder}/log.txt")

gmsh.finalize()



