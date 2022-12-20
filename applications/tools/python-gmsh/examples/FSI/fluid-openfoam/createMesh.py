import gmsh
import GmshModel

def createGmshModel(name):
    # geometrical parameters
    x0 = -3.
    x1 = -0.05
    x2 = 0.05
    x3 = 3.
    y0 = 0.
    y1 = 1.
    y2 = 4.
    z0 = 0.
    z1 = 1.
    n_z_Direction = 1

    model = gmsh.model()
    lc = 0.08 # mesh size at a point
    lc2 = 0.1
    model.geo.addPoint(x0, y0, z0, lc2, 1)
    model.geo.addPoint(x1, y0, z0, lc, 2)
    model.geo.addPoint(x1, y1, z0, lc, 3)
    model.geo.addPoint(x2, y1, z0, lc, 4)
    model.geo.addPoint(x2, y0, z0, lc, 5)
    model.geo.addPoint(x3, y0, z0, lc2, 6)
    model.geo.addPoint(x3, y2, z0, lc2, 7)
    model.geo.addPoint(x0, y2, z0, lc2, 8)

    model.geo.addLine(1, 2, 1)
    model.geo.addLine(2, 3, 2)
    model.geo.addLine(3, 4, 3)
    model.geo.addLine(4, 5, 4)
    model.geo.addLine(5, 6, 5)
    model.geo.addLine(6, 7, 6)
    model.geo.addLine(7, 8, 7)
    model.geo.addLine(8, 1, 8)
    model.geo.addCurveLoop([1, 2, 3, 4, 5, 6, 7, 8], 9)
    model.geo.addPlaneSurface([9], 1)

    dthick = [i/n_z_Direction for i in range(1,n_z_Direction+1)]
    nemEle =[1]*n_z_Direction
    #print(dthick,nemEle)
    out = model.geo.extrude([(2,1)],0,0,z1-z0,nemEle,heights=dthick,recombine=True)
    print(out)

    model.geo.addPhysicalGroup(3, [out[1][1]], 1)
    model.setPhysicalName(3,1,"Volume")

    model.geo.addPhysicalGroup(2, [1,out[0][1]], 111)
    model.setPhysicalName(2,111,"frontAndBack")

    model.geo.addPhysicalGroup(2, [out[2][1],out[6][1]], 112)
    model.setPhysicalName(2,112,"lowerWall")

    model.geo.addPhysicalGroup(2, [out[3][1],out[4][1],out[5][1]], 113)
    model.setPhysicalName(2,113,"flap")

    model.geo.addPhysicalGroup(2, [out[7][1]], 114)
    model.setPhysicalName(2,114,"outlet")

    model.geo.addPhysicalGroup(2, [out[8][1]], 115)
    model.setPhysicalName(2,115,"upperWall")

    model.geo.addPhysicalGroup(2, [out[9][1]], 116)
    model.setPhysicalName(2,116,"inlet")

    model.geo.synchronize()
    #model.mesh.setAlgorithm(2,1,8)
    #model.mesh.setRecombine(2,1)
    model.mesh.generate(model.getDimension())
    #model.mesh.setOrder(2)
    return model

boundaryConfig = {
        "frontAndBack": "empty",
        "lowerWall": "wall",
        "flap": "wall",
        "outlet": "patch",
        "upperWall": "wall",
        "inlet": "patch",
        }

gmsh.initialize()
model = createGmshModel("Flap")
gmsh.write("Flap.msh")
gmshModel = GmshModel.GmshModel(model,)
gmshModel.createPolyMesh("constant/polyMesh",boundaryConfig)
gmsh.finalize()
