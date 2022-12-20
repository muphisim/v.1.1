#class to create input file
#Author(s): Van Dung NGUYEN, 2022
#
from alive_progress import alive_bar
import Model
import Configure
import numpy as np

def createInputFile(fileName, model, configure):
    """create an input file from model and configure
        Args:
            None
        Returns:
            The file is automatically saved in the current folder
        """
    
    if not(isinstance(model,Model.Model)):
        raise TypeError("Model type is not correct!")

    if not(isinstance(configure,Configure.Configure)):
        raise TypeError("Configure type is not correct!")

    mshParts = configure.getPartionsByElementGroups()
    FEM = configure.getFEMConfigure()
    MM = configure.getMMConfigure()
    solvers = configure.getSolversConfigure()
    materials = configure.getMaterialsConfigure()
    NeumannBC = configure.getNeumannBCsConfigure()
    initialBCs = configure.getInitialBCsConfigure()
    
    # start writing file
    oxFile = open(fileName,"w")
    oxFile.write("*Node\n")
    print("write node data to file: ",fileName)
    dim = model.getDimension()
    allNodes = model.getAllNodes()
    numNodes = len(allNodes)
    with alive_bar(numNodes) as bar:
        for ii in range(numNodes):
            if dim == 3:
                oxFile.write("%d,%.16g,%.16g,%.16g\n"%(ii+1,allNodes[ii][0],allNodes[ii][1],allNodes[ii][2]))
            elif dim == 2:
                oxFile.write("%d,%.16g,%.16g\n"%(ii+1,allNodes[ii][0],allNodes[ii][1]))
            elif dim == 1:
                oxFile.write("%d,%.16g\n"%(ii+1,allNodes[ii][0]))
            else:
                raise NotImplementedError(f'Dimension {dim} is not defined!')
            bar()
    oxFile.flush()
    
    print("write element data to file: ",fileName)
    oxFile.write("*Element, type=%s\n"%(model.getElementType()))
    allElements = model.getAllElements()
    numElements = len(allElements)
    with alive_bar(numElements) as bar:
        for ii in range(numElements):
            oxFile.write("%d"%(ii+1))
            for nn in allElements[ii]:
                oxFile.write(",%d"%(nn))   
            oxFile.write("\n")
            bar()
    oxFile.flush()
    
    partitions = model.getPartitions()
    numPart =len(partitions) 
    if numPart > 0:
        print("write partition data to file: ",fileName)
        for ipar in range(numPart):
            nodesPar = partitions[ipar+1][0]
            oxFile.write("*NSET=region%d\n*CPU=%d\n"%(ipar+1,ipar+1))
            for nn in nodesPar:
                oxFile.write("%d\n"%(nn))
        for ipar in range(numPart):
            elePar = partitions[ipar+1][1]
            oxFile.write("*ELSET=region%d\n*CPU=%d\n"%(ipar+1,ipar+1))
            for ee in elePar:
                oxFile.write("%d\n"%(ee)) 
    
    #partition by physical
    if mshParts is not None:
        numPart = len(mshParts.keys())
        for ipar in range(numPart):
            if ipar+1 not in mshParts.keys():
                raise RuntimeError("wrong msh part dictionary")
            nodesPar = mshParts[ipar+1][0]
            oxFile.write("*NSET=region%d\n*CPU=%d\n"%(ipar+1,ipar+1))
            for nn in nodesPar:
                oxFile.write("%d\n"%(nn))
        for ipar in range(numPart):
            elePar = mshParts[ipar+1][1]
            oxFile.write("*ELSET=region%d\n*CPU=%d\n"%(ipar+1,ipar+1))
            for ee in elePar:
                oxFile.write("%d\n"%(ee)) 
    
    
    
    print("write FEM and MM support to file: ",fileName)
    for (key, vals) in FEM.items():
        if key.upper()  == "SUPPORT":
            if vals== "All":
                oxFile.write("*FEM Nodes\nAll\n*FEM Elements\nAll\n")
            else:
                # get all nodes and element in each item
                FEMEles= set()
                FEMNodes= set()
                for item in vals:
                    elementTags = model.getGroupOfElements(elementGroupIndex=(item[0],item[1]))
                    for iel in elementTags:
                        FEMEles.add(iel)
                        
                    nodeTags = model.getGroupOfNodes(nodeGroupIndex = (item[0],item[1]))
                    for inod in nodeTags:
                        FEMNodes.add(inod)
                oxFile.write("*FEM Nodes\n")
                for inod in FEMNodes:
                    oxFile.write("%d\n"%inod)
                oxFile.write("*FEM Elements\n")
                for iel in FEMEles:
                    oxFile.write("%d\n"%iel)
        else:
            raise NotImplementedError(f'keyword {key} is not defined!')
    oxFile.flush()
    
    
    for (key, vals) in MM.items():
        if key.upper()  == "SUPPORT":
            if vals  == "All":
                oxFile.write("*MM Nodes\nAll\n*MM Elements\nAll\n")
            else:
                MMEles= set()
                MMNodes= set()
                for item in vals:
                    elementTags = model.getGroupOfElements(elementGroupIndex=(item[0],item[1]))
                    for iel in elementTags:
                        MMEles.add(iel)
                    nodeTags = model.getGroupOfNodes(nodeGroupIndex = (item[0],item[1]))
                    for inod in nodeTags:
                        MMNodes.add(inod)
                oxFile.write("*MM Nodes\n")
                for inod in MMNodes:
                    oxFile.write("%d\n"%inod)
                oxFile.write("*MM Elements\n")
                for iel in MMEles:
                    oxFile.write("%d\n"%iel)
        else:
            raise NotImplementedError(f'keyword {key} is not defined!')
    oxFile.flush()
    
    
    print("write solver options to file: ",fileName)
    for sol in solvers:
        allBCs = {}
        keyFSI = ""
        keyExtraction = ""
        for key, val in sol.items():
            if key.upper() == "SOLVER":
                oxFile.write("*SOLVER, %s\n"%val.upper())
            elif key.upper() == "NUMBER OF STEPS":
                oxFile.write("*Number Of Steps\n%d\n"%val)    
            elif key.upper() == "SCALE FACTOR":
                oxFile.write("*Scale Factor\n%.16g\n"%val)
            elif key.upper() == "ADAPTIVE TIME STEP":
                oxFile.write("*ADAPTIVE TIME STEP\n")
                for vvi in val:
                    if val.index(vvi)==0:
                        oxFile.write(f"{vvi}")
                    else:
                        oxFile.write(f",{vvi}")
                oxFile.write("\n")    
            elif key.upper()=="MAX TIME STEP":
                oxFile.write("*Max Time Step\n%.16g\n"%val)
            elif key.upper()=="ALLOWABLE NUMBER OF FAILED STEPS":
                oxFile.write("*Allowable Number of Failed Steps\n%d\n"%val)
            elif key.upper() == "OUTPUTS":
                oxFile.write("*Outputs\n")
                for ss in val:
                    for s in ss:
                        if ss.index(s) >0:
                            oxFile.write(",")
                        oxFile.write(str(s))
                    oxFile.write("\n")
            elif key.upper() == "TIME":
                oxFile.write("*Time\n%.16g,%.16g\n"%(val[0],val[1]))
            elif key.upper() == "ABSOLUTE TOLERANCE":
                oxFile.write("*Absolute Tolerance\n%.16g\n"%(val))
            elif key.upper() == "RELATIVE TOLERANCE":
                oxFile.write("*Relative Tolerance\n%.16g\n"%(val))
            elif key.upper()=="MAX ITERATIONS":
                oxFile.write("*Max Iterations\n%d\n"%(val))
            elif key.upper() == "MAX TIME STEP":
                oxFile.write("*Max Time Step\n%d\n"%(val))
            elif key.upper() == "ERROR ANALYSIS":
                oxFile.write("*ERROR ANALYSIS\n")
                oxFile.write(f"{val[0]}")
                for iiv in range(1,len(val)):
                    oxFile.write(f",{val[iiv]}")
                oxFile.write("\n")
            elif key.upper() == "DISPBC" or key.upper() == "NODEDISPBC":
                for bcItem in val:
                    bcType = bcItem[0].upper()
                    if bcType not in allBCs.keys():
                        allBCs[bcType] = {}
                    if "DISPLACEMENT"  in bcType:
                        location = bcItem[1]
                        comVal = bcItem[2]
                        support = allBCs[bcType]
                        if key.upper() == "DISPBC":
                            nodeTags = model.getGroupOfNodes(nodeGroupIndex = (location[0],location[1]))
                        elif key.upper() == "NODEDISPBC":
                            nodeTags = [location]
                        for nn in nodeTags:
                            for cv in comVal:
                                if cv[0] not in support:
                                    support[cv[0]]={}
                                support[cv[0]][nn] = cv[1]
            elif key.upper() == "FSI":
                keyFSI = key
            elif key.upper() == "EXTRACTION":
                keyExtraction = key
    
        for key, val in allBCs.items():
            if "DISPLACEMENT"  in key:
                 oxFile.write("*BOUNDARY, TYPE=%s\n"%key.upper())
                 for com, support in val.items():
                     for node, value in support.items():
                         oxFile.write("%d,%d,%.16g\n"%(node,com,value))
                         
        if keyFSI in sol.keys():
            FSINodes= set()
            allLocations = sol[keyFSI]
            for location in allLocations:
                nodeTags = model.getGroupOfNodes(nodeGroupIndex =(location[0],location[1]))
                for inod in nodeTags:
                    FSINodes.add(inod)
            if len(FSINodes) >0:
                oxFile.write("*FSI\n")
                for inod in FSINodes:
                    oxFile.write("%d\n"%inod)
        if keyExtraction in sol.keys():
            allExtraction = sol[keyExtraction]
            for elem in allExtraction:
                if elem[0] =="NODE":
                    prefix= elem[1]
                    quantity =elem[2]
                    operation = elem[3]
                    loc = elem[4]
                    comp = elem[5]
                    oxFile.write("*EXTRACTION, NODE\n%s,%s,%s,%d\n"%(prefix,quantity,operation,comp))
                    if isinstance(loc, (list, tuple, np.ndarray)):
                        nodeTags = model.getGroupOfNodes(nodeGroupIndex =(loc[0],loc[1]))
                    else:
                        nodeTags = [loc]
                        if not(model.nodeInModel(loc)):
                            raise NotImplementedError(f'Node {loc} does not exists')
                    for inod in nodeTags:
                        oxFile.write("%d\n"%inod)
                elif elem[0] == "ELEMENT":
                    prefix= elem[1]
                    quantity =elem[2]
                    operation = elem[3]
                    loc = elem[4]
                    oxFile.write("*EXTRACTION, ELEMENT\n%s,%s,%s\n"%(prefix,quantity,operation))
                    if loc == "All":
                        oxFile.write("All\n")
                    else:
                        nodeTags = model.getGroupOfElements(elementGroupIndex =(loc[0],loc[1]))
                        for inod in nodeTags:
                            oxFile.write("%d\n"%inod)
                else:
                    raise NotImplementedError(f'keyword {elem[0]} is not defined!')
        if "PETSC_SOLVER" in sol.keys():
            oxFile.write("*PETSC_SOLVER, %s\n"%sol["PETSC_SOLVER"])
        if "PETSC_PC" in sol.keys():
            oxFile.write("*PETSC_PC, %s\n"%sol["PETSC_PC"])    
        oxFile.write("*END SOLVER\n")
        
    print("write Neumann BCs to file: ",fileName)
    for key, vals in NeumannBC.items():
        print("start writting ",key)
        if key.upper()=="GRAVITY":
           oxFile.write("*GRAVITY\n%.16g\n"%(vals)) 
        else:
            for supp in vals:
                location = supp[1]
                t0 = supp[2][0]
                tf = supp[2][1]
                values = supp[3]
                parentTags = []
                if location != "ALL" and key.upper()!="ELEMENTPRESSUREBC":
                    if location[0] ==  dim-1:
                        parentTags = model.getGroupOfFacets(facetGroupIndex=(location[0],location[1]))
                    elif location[0] == dim:
                        parentTags = model.getGroupOfElements(elementGroupIndex =(location[0],location[1]))
                    else:
                        raise NotImplementedError(f'dim {dim} cannot be used!')    
                if key.upper()=="PRESSUREBC":
                    for prEl in parentTags:
                        loc = prEl[0]
                        support = prEl[1]
                        oxFile.write("*BOUNDARY, TYPE=%s\n"%supp[0].upper())
                        oxFile.write("%.16g,%.16g,%.16g,%d\n"%(t0,tf,values,loc))
                        for jj in support:
                            oxFile.write("%d\n"%jj)
                elif key.upper()=="ELEMENTPRESSUREBC":
                    loc = location[1]
                    support = location[0]
                    oxFile.write("*BOUNDARY, TYPE=%s\n"%supp[0].upper())
                    oxFile.write("%.16g,%.16g,%.16g,%d\n"%(t0,tf,values,loc))
                    oxFile.write("%d\n"%support)
                elif key.upper()=="TRACTIONBC":                        
                    for prEl in parentTags:
                        loc = prEl[0]
                        support = prEl[1]
                        oxFile.write("*BOUNDARY, TYPE=%s\n"%supp[0].upper())
                        oxFile.write("%.16g,%.16g,%.16g,%.16g,%.16g,%d\n"%(t0,tf,values[0],values[1],values[2],loc))
                        for jj in support:
                            oxFile.write("%d\n"%jj)
                elif key.upper() == "FLUX EXTRADOF" and supp[0] == "HEAT INST":
                    for prEl in parentTags:
                        loc = prEl[0]
                        support = prEl[1]
                        oxFile.write("*FLUX EXTRADOF, TYPE=HEAT INST\n")
                        oxFile.write("%.16g,%.16g,%.16g,%.16g,%.16g,%d\n"%(t0,tf,values[0],values[1],values[2],loc))
                        for jj in support:
                            oxFile.write("%d\n"%jj)
                elif key.upper() == "FLUX EXTRADOF" and  supp[0]== "VOLUMETRIC HEAT FLUX":        
                    oxFile.write("*FLUX EXTRADOF, TYPE=VOLUMETRIC HEAT FLUX\n")
                    oxFile.write("%.16g,%.16g,%.16g\n"%(t0,tf,values))
                    for jj in parentTags:
                        oxFile.write("%d\n"%jj)   
                elif key.upper() == "FLUX EXTRADOF" and  supp[0]== "GAUSSIAN VOLUMETRIC HEAT FLUX":        
                    oxFile.write("*FLUX EXTRADOF, TYPE=GAUSSIAN VOLUMETRIC HEAT FLUX\n")
                    oxFile.write("%.16g,%.16g"%(t0,tf))
                    for vx in values:
                        oxFile.write(",%.16g"%vx)
                    oxFile.write("\n")                              
                elif key.upper() == "FLUX EXTRADOF" and  supp[0] == "CONVECTION":
                    for prEl in parentTags:
                        loc = prEl[0]
                        support = prEl[1]

                        oxFile.write("*FLUX EXTRADOF, TYPE=CONVECTION\n")
                        oxFile.write("%.16g,%.16g"%(t0,tf))
                        for cc in values:
                            oxFile.write(",%s"%str(cc))
                        oxFile.write(",%d\n"%loc)
                        for jj in support:
                            oxFile.write("%d\n"%jj)
                elif key.upper() == "FLUX EXTRADOF" and  supp[0] == "RADIATION":
                    for prEl in parentTags:
                        loc = prEl[0]
                        support = prEl[1]
                        oxFile.write("*FLUX EXTRADOF, TYPE=RADIATION\n")
                        oxFile.write("%.16g,%.16g"%(t0,tf))
                        for cc in values:
                            oxFile.write(",%s"%str(cc))
                        oxFile.write(",%d\n"%loc)
                        for jj in support:
                            oxFile.write("%d\n"%jj)
                                    
                else:
                    raise NotImplementedError(f'keyword {key.upper()} and {supp[0]} is not defined!')
        print("done writting ",key)
        
    print("write material laws to file: ",fileName)
    for mat in materials:
        for key, vals in mat.items():
            if key.upper() == "NAME":
                oxFile.write("*MATERIAL, name=%s\n"%vals)
            elif key.upper() == 'DENSITY':
                oxFile.write("*Density\n%.16g\n"%(vals))
            elif key.upper() == 'TYPE':
                oxFile.write("*%s\n"%vals)
            elif key.upper() == "PARAMETERS":
                for rowVals in vals:
                     for vv in rowVals:
                         if rowVals.index(vv) == 0:
                             oxFile.write('%s'%str(vv))
                         else:
                             oxFile.write(',%s'%str(vv))
                     oxFile.write("\n")
            elif key.upper() == "EXTRADOF":
                oxFile.write("*ExtraDof\n")
                oxFile.write("*%s\n"%(vals[0]))
                for ii in range(len(vals[1])):
                    if ii >0:
                        oxFile.write(",")
                    oxFile.write("%s"%str(vals[1][ii]))
                oxFile.write("\n")
            elif key.upper() == "LIST OF ELEMENTS":
                if vals== "All":
                    oxFile.write("*List of Elements\nAll\n")
                else:
                    allEles = set()
                    for item in vals:
                        elementTags =model.getGroupOfElements(elementGroupIndex = (item[0],item[1]))
                        for iel in elementTags:
                            allEles.add(iel)
                    oxFile.write("*List of Elements\n")
                    for iel in allEles:
                        oxFile.write("%d\n"%iel)
            else:
                raise NotImplementedError(f'keyword {key} is not defined!')
    for key, vals in initialBCs.items():
        if key.upper() == "INITIAL CONDITIONS EXTRADOF":
            oxFile.write("*INITIAL CONDITIONS EXTRADOF\n")
            for ii in vals:
                if vals.index(ii) >0:
                    oxFile.write(",")
                oxFile.write(str(ii))
            oxFile.write("\n")    
        
    oxFile.write("*END")
    print("done writing input file: ",fileName)
    oxFile.close()
    
