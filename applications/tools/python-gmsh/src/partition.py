#!/usr/bin/python3

import sys
import gmsh
import os

a= sys.argv
b = " ".join(a)
print("run command with python: "+b)
if len(a) ==4:
  inputFile = a[1]
  gmshFileName = inputFile.split('.')[0]+'-temp.msh'
  inputFileNew = a[2]
  nparts = int(a[3])
  print(f"create new input file with {nparts} parts")
else:
  raise SyntaxError("use: partition.y inputFile outputFile nparts")
  sys.exit()


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

gmshFile = open(gmshFileName, "w")
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
    if nodeListCreated and eleListCreated:
        break
    else:
        lineIndex = lineIndex+1
        if lineIndex >= len(Lines):
            break


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
gmshFile.write("$EndNodes\n$Elements\n%d\n"%(numElements));
for iele in range(len(ElementList)):
    llNoSpace=ElementList[iele].strip()
    linspl = llNoSpace.split(',')
    phyGrp = 1
    eleID = int(linspl[0])
    gmshFile.write('%s %s 2 %d %d'%(linspl[0],AbaqusToGmsh[meshType],phyGrp,phyGrp))
    linspl.pop(0)
    convertNodeNumbering(AbaqusToGmsh[meshType],linspl)
    for jj in linspl:
        gmshFile.write(" %s"%jj)
    gmshFile.write("\n")    
gmshFile.write("$EndElements");
gmshFile.close()

gmsh.initialize()
gmsh.open(gmshFileName)

partitions = {}
if nparts >1:
    gmsh.model.mesh.partition(nparts)
    entities = gmsh.model.getEntities()
    for e in entities:
        if e[0] == gmsh.model.getDimension() :
            procs = gmsh.model.getPartitions(e[0], e[1])
            elementTypes, elementTags, nodeTags =  gmsh.model.mesh.getElements(e[0], e[1])
            for par in procs:
                if par not in partitions.keys():
                    partitions[par] = (set(),set())
                nodeSet = partitions[par][0]
                for nno in nodeTags:
                    for nn in nno:
                        nodeSet.add(nn)
                elSet = partitions[par][1]
                for eleo in elementTags:
                    for ee in eleo:
                        elSet.add(ee)

newFile = open(inputFileNew, "w")
for ll in Lines[:lineIndex+1]:
    newFile.write(ll)
if nparts > 1:
    for ipar in range(nparts):
        nodesPar = partitions[ipar+1][0]
        newFile.write("*NSET=region%d\n*CPU=%d\n"%(ipar+1,ipar+1))
        for nn in nodesPar:
            newFile.write("%d\n"%(nn))
    for ipar in range(nparts):
        elePar = partitions[ipar+1][1]
        newFile.write("*ELSET=region%d\n*CPU=%d\n"%(ipar+1,ipar+1))
        for ee in elePar:
            newFile.write("%d\n"%(ee))                 

for ll in Lines[lineIndex+1:]:
    newFile.write(ll)
    
newFile.close()
gmsh.finalize()

os.system("rm -rf "+gmshFileName)
