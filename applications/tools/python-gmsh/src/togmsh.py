#!/usr/bin/python

import sys
import gmsh
import os

a= sys.argv
b = " ".join(a)
print("run command with python: "+b)
if len(a) >1:
  inputFile = a[1]
  if (len(a)==2):
    outputFile = inputFile.split('.')[0]+'-temp.msh'
  else:
    outputFile = a[2]

  print('Input file: %s\nOutput file: %s'%(inputFile,outputFile))
else:
  sys.exit("inputFile, outputFile must be provided")


AbaqusToGmsh={'None':'0',
        'C3D4':'4',
        'C3D4H':'4',
        'CPS3':'2',
        'CPE3':'2',
        'CPS4':'3',
        'CPE4':'3',
        'CPS4R':'3',
        'CPS6':'9',
        'C3D8':'5',
        'C3D8R':'5',
        'C3D10':'11',
        'B21':'1',
        'B22':'8',
        'B31':'1',
        'B32':'8',
        }
ElementDimension ={'None':'-1',
        'C3D4':'3',
        'C3D4H':'3',
        'CPS3':'2',
        'CPE3':'2',
        'CPS4':'2',
        'CPE4':'2',
        'CPS4R':'2',
        'CPS6':'2',
        'C3D8':'3',
        'C3D8R':'3',
        'C3D10':'3',
        'B21':'1',
        'B22':'1',
        'B31':'1',
        'B32':'1',
        }

def convertNodeNumbering(gmshType, nodeList):
        """Node numberfing in gmsh and Abaqus are diffrent,
        this function is used to convert from gmsh node-ordering to Abaqus node ordering

        """
        if gmshType == "11":
            # 10 nodes tetrahedron
            tmp = nodeList[9]
            nodeList[9] = nodeList[8]
            nodeList[8] = tmp
        elif gmshType == "8":
            # 10 nodes tetrahedron
            tmp = nodeList[2]
            nodeList[2] = nodeList[1]
            nodeList[1] = tmp

gmshFile = open(outputFile, "w")
gmshFile.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")


ff = open(inputFile, "r",encoding='latin1')
Lines = ff.readlines()
ff.close()

lineIndex=0
nodeListCreated=False
NodeList = []
eleListCreated=False
ElementList = []
meshType='None'
NodeGroups = []
ElementGroups = {}
NodeGroups = {}

printTimeGR={}
invarsTimeGr = {}
displacementBCTimeGr = {}

while True:
    line = Lines[lineIndex]
    line=line.strip()
    if line[:5].upper()=="*NODE" and not(nodeListCreated):
        nodeListCreated = True
        while True:
            lineIndex = lineIndex+1
            if lineIndex >= len(Lines):
                break
            line = Lines[lineIndex]
            line=line.strip()
            if line[0] == "*":
                lineIndex = lineIndex-1
                break
            else:
                NodeList.append(line)

    elif line[:8].upper()=="*ELEMENT" and not(eleListCreated):
        aa = line.split("=")
        meshType=aa[-1]
        eleListCreated = True
        while True:
            lineIndex = lineIndex+1
            if lineIndex >= len(Lines):
                break
            line = Lines[lineIndex]
            line=line.strip()
            if line[0] == "*":
                lineIndex = lineIndex-1
                break
            else:
                ElementList.append(line)

    elif line[:6].upper() == "*ELSET":
        aa = line.split("=",1)
        physicalName=aa[-1]
        ElementGroups[physicalName] = []
        lineNoSpace = line.replace(" ","")
        lineNoSpace = lineNoSpace.strip()
        lineSplitted = lineNoSpace.split(",")
        while True:
            lineIndex = lineIndex+1
            if lineIndex >= len(Lines):
                break
            line = Lines[lineIndex]
            line=line.strip()
            if line[0] == "*":
                lineIndex = lineIndex-1
                break
            else:
                line = line.replace(","," ")
                line = line.strip()
                lineSp = line.split()
                if lineSplitted[-1].upper() == "GENERATE":
                    startv = int(lineSp[0].strip())
                    endv =  int(lineSp[1].strip())
                    incre = int(lineSp[2].strip())
                    for eleNum in range(startv,endv+1,incre):
                        ElementGroups[physicalName].append(eleNum)
                else: 
                    for ee in lineSp:
                        eleNum = int(ee)
                        ElementGroups[physicalName].append(eleNum)
                        
    elif line[:5].upper() == "*NSET":
        aa = line.split("=",1)
        physicalName=aa[-1]
        NodeGroups[physicalName] = []
        lineNoSpace = line.replace(" ","")
        lineNoSpace = lineNoSpace.strip()
        lineSplitted = lineNoSpace.split(",")
        
        while True:
            lineIndex = lineIndex+1
            if lineIndex >= len(Lines):
                break
            line = Lines[lineIndex]
            line=line.strip()
            if line[0] == "*":
                lineIndex = lineIndex-1
                break
            else:
                line = line.replace(","," ")
                line = line.strip()
                lineSp = line.split()
                if lineSplitted[-1].upper() == "GENERATE":
                    startv = int(lineSp[0].strip())
                    endv =  int(lineSp[1].strip())
                    incre = int(lineSp[2].strip())
                    for eleNum in range(startv,endv+1,incre):
                        NodeGroups[physicalName].append(eleNum)
                else: 
                    for ee in lineSp:
                        eleNum = int(ee)
                        NodeGroups[physicalName].append(eleNum)
    elif line[:28].upper() == "*BOUNDARY, TYPE=DISPLACEMENT":
        prefix = line[16:]
        while True:
            lineIndex = lineIndex+1
            if lineIndex >= len(Lines):
                break
            line = Lines[lineIndex]
            if line[0] == "*":
                lineIndex = lineIndex-1
                break
            else:
                line = line.replace(","," ")
                line = line.strip()
                lineSp = line.split()
                nodeNum = int(lineSp[0])
                comp = int(lineSp[1])
                value = float(lineSp[2])
                
                found = False
                for prGrName, vv in displacementBCTimeGr.items():
                    if vv[0] == comp and  abs(vv[1]-value) < 1e-8:
                        NodeGroups[prGrName].append(nodeNum)
                        found = True
                        break
                if not(found):
                    name = prefix+f"-{comp}-{value}"
                    displacementBCTimeGr[name] = [comp,value]    
                    NodeGroups[name]=[nodeNum]
    elif line[:41].upper() == "*FLUX EXTRADOF, TYPE=VOLUMETRIC HEAT FLUX":
        eleList = []
        lineIndex = lineIndex+1
        line = Lines[lineIndex]
        line = line.replace(","," ")
        line = line.strip()
        lineSp = line.split()
        startTime = float(lineSp[0])
        endTime = float(lineSp[1])
        lsPower = float(lineSp[2])
        while True:
            lineIndex = lineIndex+1
            if lineIndex >= len(Lines):
                break
            line = Lines[lineIndex]
            if line[0] == "*":
                lineIndex = lineIndex-1
                break
            else:
                line = line.replace(","," ")
                line = line.strip()
                lineSp = line.split()
                eleList.append(int(lineSp[0]))
        
        found = False
        for prGrName, vv in printTimeGR.items():
            if abs(vv[0]-startTime)<1e-8 and  abs(vv[1]-endTime) < 1e-8:
                found =True
                for ee in eleList:
                    ElementGroups[prGrName].append(ee)
                break
        if not(found):
            name = f"LaserPower-{startTime:g}-{endTime:g}-{lsPower:g}"
            printTimeGR[name] = [startTime,endTime]
            ElementGroups[name] = eleList
                  
    elif line[:12].upper() == "*INITINTVARS":        
        while True:
            lineIndex = lineIndex+1
            if lineIndex >= len(Lines):
                break
            line = Lines[lineIndex]
            if line[0] == "*":
                lineIndex = lineIndex-1
                break
            else:
                line = line.replace(","," ")
                line = line.strip()
                lineSp = line.split()
                eleNum = int(lineSp[0])
                comp = int(lineSp[1])
                value = float(lineSp[2])
                
                found = False
                for prGrName, vv in invarsTimeGr.items():
                    if vv[0] == comp and  abs(vv[1]-value) < 1e-8:
                        ElementGroups[prGrName].append(eleNum)
                        found = True
                        break
                if not(found):
                    name = f"INITINTVARS-{comp}-{value}"
                    invarsTimeGr[name] = [comp,value]    
                    ElementGroups[name] = [eleNum]
                    
                             
    lineIndex = lineIndex+1
    if lineIndex >= len(Lines):
        break

#print(NodeList)
#print(ElementList)
#print(meshType)
#print(ElementGroups.keys())
#print(NodeGroups)
#NodeGroups.clear()

numNodesInNodeGroups=0
numElementsInElementGroups =0
physicalGroupIdex={}
if len(ElementGroups)+len(NodeGroups) >0:
    gmshFile.write(f"$PhysicalNames\n{len(ElementGroups)+len(NodeGroups)}\n")
    for ii in ElementGroups.keys():
        oldLen = len(physicalGroupIdex)
        physicalGroupIdex[ii] = oldLen+1
        eleDim = ElementDimension[meshType]
        gmshFile.write(f"{eleDim} {oldLen+1} \"{ii}\"\n")
        numElementsInElementGroups += len(ElementGroups[ii])

    for ii in NodeGroups.keys():
        oldLen = len(physicalGroupIdex)
        physicalGroupIdex[ii] = oldLen+1
        numNodesInNodeGroups += len(NodeGroups[ii])
        gmshFile.write(f"0 {oldLen+1} \"{ii}\"\n")
    gmshFile.write("$EndPhysicalNames\n")

numNodes = len(NodeList)
numElements=len(ElementList)

print('Number of nodes: %d\nNumber of elements: %d\nElememt type: %s'%(numNodes,numElements,meshType))
#
gmshFile.write("$Nodes\n");
gmshFile.write("%d\n"%numNodes)
#
for inode in range(numNodes):
    llNoSpace=NodeList[inode].strip()
    linspl = llNoSpace.split(',')
    if (len(linspl) < 4):
        while len(linspl)!=4:
            linspl.append('0.0')
    #print(linspl)
    for st in linspl:
        gmshFile.write("%s "%st)
    gmshFile.write("\n")

gmshFile.write("$EndNodes\n$Elements\n%d\n"%(numElements+numElementsInElementGroups+numNodesInNodeGroups));

elementsDict = {}
for iele in range(len(ElementList)):
    llNoSpace=ElementList[iele].strip()
    linspl = llNoSpace.split(',')
    phyGrp = len(physicalGroupIdex)+1
    eleID = int(linspl[0])
    gmshFile.write('%s %s 2 %d %d'%(linspl[0],AbaqusToGmsh[meshType],phyGrp,phyGrp))
    linspl.pop(0)
    convertNodeNumbering(AbaqusToGmsh[meshType],linspl)
    for jj in linspl:
        gmshFile.write(" %s"%jj)
    gmshFile.write("\n")
    elementsDict[eleID] = linspl
    
eleStart=max(elementsDict.keys())+1
for jj, vals in ElementGroups.items():
    phyGrp = physicalGroupIdex[jj]
    for ele in vals:
        linspl = elementsDict[ele]
        gmshFile.write('%d %s 2 %d %d'%(eleStart,AbaqusToGmsh[meshType],phyGrp,phyGrp))
        for kk in linspl:
            gmshFile.write(" %s"%kk)
        gmshFile.write("\n")
        eleStart = eleStart+1
        
for jj, val in NodeGroups.items():
    phyGrp = physicalGroupIdex[jj]
    for nn in val:
        gmshFile.write('%s 15 2 %d %d %d\n'%(eleStart,phyGrp,nn,nn))
        eleStart = eleStart+1
gmshFile.write("$EndElements");

gmshFile.close()

gmsh.initialize()
gmsh.open(outputFile)
gmsh.fltk.run()
os.system("rm -rf "+outputFile)
gmsh.finalize()
