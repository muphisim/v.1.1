#class to define the mesh model with INP input file of MuPhiSim
#Author(s): Van Dung NGUYEN, 2022
#

from alive_progress import alive_bar
import os
import shutil
import numpy as np
import Model
import GmshModel

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

class InpModel (Model.Model):
    
    def __init__(self, inpInputFileName):
        """ a model from inp input file
        Args:
            inpInputFileName: name of input file
        """
        
        self.dimension=-1
        self.elementType=""
        self.nodes=[]
        self.elements=[]
        
        self.partitions={}
        self.nodeMap={}
        self.elementMap={}
        self.edges=[]
        self.groupOfElements={}
        self.groupOfNodes={}
        self.groupOfFacets={}
        
        with open(inpInputFileName, 'r') as fInput:
            # do things with your file

            Lines = fInput.readlines()
            lineIndex=0
            nodeListCreated=False
            eleListCreated=False
            ElementGroups = {}
            NodeGroups = {}
            displacementBCTimeGr = {}
            print(f"start reading file {inpInputFileName}")
            with alive_bar(len(Lines)) as bar:
                while True:
                    line = Lines[lineIndex]
                    line=line.strip()
                    if line[:5].upper()=="*NODE" and not(nodeListCreated):
                        nodeListCreated = True
                        allNodes=[]
                        while True:
                            lineIndex = lineIndex+1
                            bar()
                            if lineIndex >= len(Lines):
                                break
                            line = Lines[lineIndex]
                            line=line.strip()
                            if line[0] == "*":
                                lineIndex = lineIndex-1
                                break
                            else:
                                lineSplitted=line.split(",")
                                num = int(lineSplitted[0])
                                self.nodeMap[num]=num
                                self.dimension=len(lineSplitted[1:])
                                if self.dimension==2:
                                    allNodes.append([float(lineSplitted[1]),float(lineSplitted[2]),0])
                                elif self.dimension==3:
                                    allNodes.append([float(lineSplitted[1]),float(lineSplitted[2]),float(lineSplitted[3])])
                                else:
                                    raise RuntimeError(f"dimension {self.dimension} is not implemented")
                        self.nodes=np.array(allNodes,dtype=float)

                    elif line[:8].upper()=="*ELEMENT" and not(eleListCreated):
                        aa = line.split("=")
                        self.elementType=aa[-1].strip()
                        eleListCreated = True
                        ElementList=[]
                        while True:
                            lineIndex = lineIndex+1
                            bar()
                            if lineIndex >= len(Lines):
                                break
                            line = Lines[lineIndex]
                            line=line.strip()
                            if line[0] == "*":
                                lineIndex = lineIndex-1
                                break
                            else:
                                lineSplitted=line.split(",")
                                num=int(lineSplitted[0])
                                self.elementMap[num]=num
                                ElementList.append([int(vv) for vv in lineSplitted[1:]])
                        self.elements=np.array(ElementList,dtype=int)
                        
                    elif line[:6].upper() == "*ELSET":
                        aa = line.split("=",1)
                        physicalName=aa[-1]
                        ElementGroups[physicalName] = []
                        while True:
                            lineIndex = lineIndex+1
                            bar()
                            if lineIndex >= len(Lines):
                                break
                            line = Lines[lineIndex]
                            line=line.strip()
                            if (line[0] == "*") and (line[:4] != "*CPU"):
                                lineIndex = lineIndex-1
                                break
                            elif line[:4] == "*CPU":
                                print(line)
                            else:
                                line = line.strip()
                                lineSp = line.split()
                                for ee in lineSp:
                                    eleNum = int(ee)
                                    ElementGroups[physicalName].append(eleNum)
                  
                    elif line[:5].upper() == "*NSET":
                        aa = line.split("=",1)
                        physicalName=aa[-1]
                        NodeGroups[physicalName] = []
                        while True:
                            lineIndex = lineIndex+1
                            bar()
                            if lineIndex >= len(Lines):
                                break
                            line = Lines[lineIndex]
                            line=line.strip()
                            if (line[0] == "*") and (line[:4] != "*CPU"):
                                lineIndex = lineIndex-1
                                break
                            elif line[:4] == "*CPU":
                                print(line)
                            else:
                                line = line.strip()
                                lineSp = line.split()
                                for ee in lineSp:
                                    eleNum = int(ee)
                                    NodeGroups[physicalName].append(eleNum)
                    elif line[:28].upper() == "*BOUNDARY, TYPE=DISPLACEMENT":
                        prefix = line[16:]
                        while True:
                            lineIndex = lineIndex+1
                            bar()
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
                    elif line[0]=="*":
                        print(f"keyword {line} is not implemented yet") 
                                                
                    lineIndex = lineIndex+1
                    bar()
                    if lineIndex >= len(Lines):
                        break
                
            print(f"done reading file {inpInputFileName}")
            index=1
            for key, val in ElementGroups.items():
                self.groupOfElements[(self.dimension,index)]=[v for v in val]
                index = index+1
            
            index =1   
            for key, val in NodeGroups.items():
                self.groupOfNodes[(0,index)]=[v for v in val]
                index = index+1
                
    def writeToGmshForView(self, gmshFileName):
        gmshFile = open(gmshFileName, "w")
        gmshFile.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        
        numNodes = len(self.nodes)
        numElements=len(self.elements)
        numElementsGroups = len(self.groupOfElements)
        numNodesGroups = len(self.groupOfNodes)
        
        numNodesInNodeGroups=0
        numElementsInElementGroups =0
        if numElementsGroups+numNodesGroups >0:
            gmshFile.write(f"$PhysicalNames\n{numElementsGroups+numNodesGroups}\n")
            for kk, val in self.groupOfElements.items():
                eleDim = ElementDimension[self.elementType]
                gmshFile.write(f"{eleDim} {kk[1]} \"{kk[1]}\"\n")
                numElementsInElementGroups += len(val)

            for kk, val in self.groupOfNodes.items():
                numNodesInNodeGroups += len(val)
                gmshFile.write(f"0 {kk[1]} \"{kk[1]}\"\n")
            gmshFile.write("$EndPhysicalNames\n")
            
        
        gmshFile.write("$Nodes\n");
        gmshFile.write("%d\n"%numNodes)
        #
        for inode in range(numNodes):
            gmshFile.write("%d"%(inode+1))
            for st in self.nodes[inode]:
                gmshFile.write(" %g"%st)
            gmshFile.write("\n")
        gmshFile.write("$EndNodes\n$Elements\n%d\n"%(numElements+numElementsInElementGroups+numNodesInNodeGroups));
        
        elementsDict = {}
        for iele in range(len(self.elements)):
            phyGrp = len(self.groupOfElements)+10
            eleID = iele+1
            gmshFile.write('%s %s 2 %d %d'%(eleID,AbaqusToGmsh[self.elementType],phyGrp,phyGrp))
            linspl = [vv for vv in self.elements[iele]]
            convertNodeNumbering(AbaqusToGmsh[self.elementType],linspl)
            for jj in linspl:
                gmshFile.write(" %d"%jj)
            gmshFile.write("\n")
            elementsDict[eleID] = linspl
            
        eleStart=max(elementsDict.keys())+1
        for jj, vals in self.groupOfElements.items():
            for ele in vals:
                linspl = elementsDict[ele]
                gmshFile.write('%d %s 2 %d %d'%(eleStart,AbaqusToGmsh[self.elementType],jj[1],jj[1]))
                for kk in linspl:
                    gmshFile.write(" %d"%kk)
                gmshFile.write("\n")
                eleStart = eleStart+1
                
        for jj, val in self.groupOfNodes.items():
            for nn in val:
                gmshFile.write('%s 15 2 %d %d %d\n'%(eleStart,jj[1],nn,nn))
                eleStart = eleStart+1
        gmshFile.write("$EndElements");

        gmshFile.close()

              
    def getPartitions(self):
        """return all partitions
        Args:
            None
        Returns:
            {1: (nodeSet, elementSet), ...}
        """
        return self.partitions
    def nodeInModel(self, nodeId):
        """return a boolean
        Args:
            nodeId: Node Id
        Returns:
            True if node is in the model, false otherwise
        """
        return nodeId in self.nodeMap

    def elementInModel(self, elementId):
        """return a boolean
        Args:
            elementId: element Id
        Returns:
            True if element is in the model, false otherwise
        """
        return elementId in self.elementMap

    def getDimension(self):
        """return a dimension of the problem
        Args:
            None
            Returns:
                The dimension of the problem (2 or 3)
        """
        return self.dimension
    def getAllNodes(self):
        """return a list of all nodes in the mesh
        Args:
            None
            Returns:
                The output is an array of node coordinates
                    [[node1-x node1-y node1-z],
                     [node2-x node2-y node2-z],
                     ...
                     ]
        The row index corresponds to the node identifier-1
        """
        return self.nodes

    def getElementType(self):
        """return the element type of the mesh, only one element type is present
        Args:
            None
        Returns:
            The element type by string
        """
        return self.elementType

    def getAllElements(self):
        """return a list of all elements in the mesh
        Args:
            None
        Returns:
            The output is an array of nodes in this element
                [[node1-ele1, node2-ele1, ...],
                [node1-ele2, node2-ele2, ...],
                ...
                ]
            The row index corresponds to the element identifier-1
        """
        return self.elements
        
    def getAllEdges(self):
        """return a list of all edges of elements in the mesh
        Args:
            None
        Returns:
            The output is an array of nodes in edges
                [[node1-ele1, node2-ele1, ...],
                [node1-ele2, node2-ele2, ...],
                ...
                ]
        """
        if self.dimension == 1:
            return self.elements
        else:
            return self.edges
        
    def getGroupOfElements(self, elementGroupIndex):
        """return a list of elements in a group identified elementGroupIndex
        Args:
            elementGroupIndex: unique pair (dim, groupId) to identify the element group
        Returns:
            List of element identifiers.
        """
        if elementGroupIndex in self.groupOfElements.keys():
            return self.groupOfElements[elementGroupIndex]
        else:
            print(f"group of element {elementGroupIndex} does not exist")
            return [] # return an empty list if no group found

    def getGroupOfNodes(self, nodeGroupIndex):
        """return a list of nodes in a group identified nodeGroupIndex
        Args:
            nodeGroupIndex: unique number to identify the group of nodes
        Returns:
            List of node identifiers.
        """
        if nodeGroupIndex in self.groupOfNodes.keys():
            return self.groupOfNodes[nodeGroupIndex]
        else:
            print(f"group of nodes {nodeGroupIndex} does not exist")
            return [] # return an empty list if no group found
    
    def getGroupOfNodesIndexes(self):
        """return a list of nodeGroupIndex
        Args:
            Node
        Returns:
            List of nodeGroupIndex.
        """
        return list(self.groupOfNodes.keys())
    
    def getGroupOfElementsIndexes(self):
        """return a list of elementGroupIndex
        Args:
            Node
        Returns:
            List of elementGroupIndex.
        """
        return list(self.groupOfElements.keys())

    def getGroupOfFacets(self, facetGroupIndex):
        """return a list of facet in a group identified nodeGroupIndex
        Args:
            facetGroupIndex: unique number to identify the group of facets
        Returns:
            List of facet position and a list of element identifiers.
                [[faceposition, [ele1, ele2, ....]],
                 [faceposition, [ele1, ele2, ....]]
                ]
        """
        if facetGroupIndex in self.groupOfFacets.keys():
            return self.groupOfFacets[facetGroupIndex]
        else:
            print(f"group of facets {facetGroupIndex} does not exist")
            return [] # return an empty list if no group found
            
    def centroid_refine(self, boundaryOnly=False, elementGroups=None, level=1):
        """refine mesh form centroid
        Args:
            elementGroups: None if refinement is applied on all element; if not, the elements in this groups is refined
            boundaryOnly: True if only element close to the boundary is refined
        Returns:
            None
        """
        if self.dimension != 2:
            raise RuntimeError("this method is implemented for 2d cases only and linear triangular elements")
        
        print("start mesh refining")
        
        refineElements=[]
        if elementGroups is None:
            refineElements = [i+1 for i in range(len(self.elements))]
        else:
            refineElements = [ee for ee in elementGroups]
            
        if boundaryOnly:
            raise NotImplementedError("boundaryOnly is not implemented")
        
        numEleOld = len(self.elements)
        oldElements = np.array(self.elements)
        newElementsMap = {}
        newElements = []
        allNodes = self.nodes.tolist()
        numNode = len(self.nodes) 
        for i in range(numEleOld):
            nodeLists=self.elements[i]
            if i+1 in refineElements:
                coords = np.array(self.nodes[nodeLists-1])
                aveNew =  coords.mean(axis=0)
                allNodes.append(aveNew.tolist())
                newElements.append([nodeLists[0],nodeLists[1],len(allNodes)])
                newElements.append([nodeLists[1],nodeLists[2],len(allNodes)])
                newElements.append([nodeLists[2],nodeLists[0],len(allNodes)])
                newElementsMap[i+1]=[len(newElements),len(newElements)-1,len(newElements)-2]
            else:
                newElements.append([jj for jj in nodeLists])
                newElementsMap[i+1]=[len(newElements)]
             
        
        for key in self.groupOfElements.keys():
            newEles = []
            newNodes = set()
            for v in self.groupOfElements[key]:
                for c in newElementsMap[v]:
                    newEles.append(c)
                    if key in self.groupOfNodes.keys():
                        for nn in newElements[c-1]:
                            newNodes.add(nn)
            self.groupOfElements[key] = np.array(newEles)
            if key in self.groupOfNodes.keys():
                self.groupOfNodes[key] = np.array(list(newNodes))
        
            
        self.nodes = np.array(allNodes)
        self.elements = np.array(newElements,dtype=int)
     
        print("done mesh refining")

