#class to define the mesh model with GMSH
#Author(s): Van Dung NGUYEN, 2022
#

from alive_progress import alive_bar
import os
import shutil
import numpy as np
import Model

class GmshModel (Model.Model):
    @staticmethod
    def getType(gmshType):
        """Return the Abaqus element type from a gmsh element type
        """
        if gmshType == 4:
            # 4-node tetrahedron
            return "C3D4"
        if gmshType == 11:
            # 10-node tetrahedron
            return "C3D10"
        elif gmshType == 2:
            # 3 node triangle
            return "CPS3"
        elif gmshType == 9:
            # 6 node triangle
            return "CPS6"
        elif gmshType == 3:
            # 4-node quad
            return "CPE4"
        elif gmshType == 6:
            #6-node linear triangular prism
            return "C3D6"
        elif gmshType == 5:
            #8-node linear brick
            return "C3D8"
        elif gmshType == 1:
            return "B21"
        else:
            raise NotImplementedError(f"Element type {gmshType} has not been defined in getType")

    @staticmethod
    def convertNodeNumbering(gmshType, nodeList):
        """Node numberfing in gmsh and Abaqus are diffrent,
        this function is used to convert from gmsh node-ordering to Abaqus node ordering

        """
        if gmshType == 11:
            # 10 nodes tetrahedron
            tmp = nodeList[9]
            nodeList[9] = nodeList[8]
            nodeList[8] = tmp
    
    @staticmethod        
    def getEdges(elementType, nodeList):
        """ return sub-elements in an element
        Args:
            elementType: element type identified by Gmsh
            nodeList: list of nodes of the element
        Returns:
            List of lines
        """
        allSublists = []
        if elementType == "C3D4":
            # 4-node tetrahedron
            allLines=[[nodeList[0],nodeList[1]],
                      [nodeList[1],nodeList[2]],
                      [nodeList[2],nodeList[0]],
                      [nodeList[0],nodeList[3]],
                      [nodeList[1],nodeList[3]],
                      [nodeList[2],nodeList[3]]]
        elif elementType == "C3D10":
            # 10-node tetrahedron
            allLines=[[nodeList[0],nodeList[1],nodeList[4]],
                      [nodeList[1],nodeList[2],nodeList[5]],
                      [nodeList[2],nodeList[0],nodeList[6]],
                      [nodeList[0],nodeList[3],nodeList[7]],
                      [nodeList[1],nodeList[3],nodeList[8]],
                      [nodeList[2],nodeList[3],nodeList[9]]]
        elif elementType == "CPS3":
            # 3-node triangle
            allLines=[[nodeList[0],nodeList[1]],
                     [nodeList[1],nodeList[2]],
                     [nodeList[2],nodeList[0]]]
        elif elementType == "CPS6":
            # 6-node triangle
            allLines=[[nodeList[0],nodeList[1],nodeList[3]],
                     [nodeList[1],nodeList[2],nodeList[4]],
                     [nodeList[2],nodeList[0],nodeList[5]]]
        elif elementType == "CPE4":
            # 4-node quad
            allLines=[[nodeList[0],nodeList[1]],
                     [nodeList[1],nodeList[2]],
                     [nodeList[2],nodeList[3]],
                     [nodeList[3],nodeList[0]]]
        elif elementType == "C3D6":
            #6-node linear triangular prism
            allLines=[[nodeList[0],nodeList[1]],
                     [nodeList[1],nodeList[2]],
                     [nodeList[2],nodeList[0]],
                     [nodeList[3],nodeList[4]],
                     [nodeList[4],nodeList[5]],
                     [nodeList[5],nodeList[3]],
                     [nodeList[0],nodeList[3]],
                     [nodeList[1],nodeList[4]],
                     [nodeList[2],nodeList[5]],]
        elif elementType == "C3D8":
            #8-node linear brick
            allLines=[[nodeList[0],nodeList[1]],
                     [nodeList[1],nodeList[2]],
                     [nodeList[2],nodeList[3]],
                     [nodeList[3],nodeList[0]],
                     [nodeList[4],nodeList[5]],
                     [nodeList[5],nodeList[6]],
                     [nodeList[6],nodeList[7]],
                     [nodeList[7],nodeList[4]],
                     [nodeList[0],nodeList[4]],
                     [nodeList[1],nodeList[5]],
                     [nodeList[2],nodeList[6]],
                     [nodeList[3],nodeList[7]],]
        else:
            raise NotImplementedError(f"Element type {elementType} has not been defined in getLines")
        return allLines

    @staticmethod
    def getFacets(elementType, nodeList):
        """ return sub-elements in an element
        Args:
            elementType: a string of element (Abaqus keyword)
            nodeList: list of nodes of the element
        Returns:
            List of facets as Abaqus
        """
        allSublists = []
        if elementType == "C3D4":
            # 4-node tetrahedron
            allSublists=[[nodeList[2],nodeList[1],nodeList[0]],
                         [nodeList[3],nodeList[0],nodeList[1]],
                         [nodeList[3],nodeList[2],nodeList[1]],
                         [nodeList[3],nodeList[2],nodeList[0]]]
        elif elementType == "C3D10":
            # 10-node tetrahedron
            allSublists=[[nodeList[2],nodeList[1],nodeList[0],nodeList[4],nodeList[5],nodeList[6]],
                         [nodeList[3],nodeList[0],nodeList[1],nodeList[7],nodeList[4],nodeList[8]],
                         [nodeList[3],nodeList[2],nodeList[1],nodeList[9],nodeList[5],nodeList[8]],
                         [nodeList[3],nodeList[2],nodeList[0],nodeList[9],nodeList[6],nodeList[7],]]
        elif elementType == "CPS3":
            # 3-node triangle
            allSublists=[[nodeList[0],nodeList[1]],
                         [nodeList[1],nodeList[2]],
                         [nodeList[2],nodeList[0]]]
        elif elementType == "CPS6":
            # 6-node triangle
            allSublists=[[nodeList[0],nodeList[1],nodeList[3]],
                         [nodeList[1],nodeList[2],nodeList[4]],
                         [nodeList[2],nodeList[0],nodeList[5]]]
        elif elementType == "CPE4":
            # 4-node quad
            allSublists=[[nodeList[0],nodeList[1]],
                         [nodeList[1],nodeList[2]],
                         [nodeList[2],nodeList[3]],
                         [nodeList[3],nodeList[0]]]
        elif elementType == "C3D6":
            #6-node linear triangular prism
            allSublists=[[nodeList[0],nodeList[2],nodeList[1]],
                         [nodeList[3],nodeList[4],nodeList[5]],
                         [nodeList[0],nodeList[1],nodeList[4],nodeList[3]],
                         [nodeList[1],nodeList[2],nodeList[5],nodeList[4]],
                         [nodeList[0],nodeList[3],nodeList[5],nodeList[2]]]
        elif elementType == "C3D8":
            #8-node linear brick
            allSublists=[[nodeList[0],nodeList[3],nodeList[2],nodeList[1]],
                         [nodeList[4],nodeList[5],nodeList[6],nodeList[7]],
                         [nodeList[0],nodeList[1],nodeList[5],nodeList[4]],
                         [nodeList[1],nodeList[2],nodeList[6],nodeList[5]],
                         [nodeList[3],nodeList[7],nodeList[6],nodeList[2]],
                         [nodeList[0],nodeList[4],nodeList[7],nodeList[3]]]
        else:
            raise NotImplementedError(f"Element type {elementType} has not been defined in getAllSubElementNodes")
        return allSublists

    def __init__(self, gmshModel, nparts=1, createFacets=False, createEdges=False):
        """ a model from gmsh
        Args:
            gmshModel: a gmsh.model object
        """
        print("start loading gmsh model")
        self.gmshModel = gmshModel
        self.dimension =gmshModel.getDimension() # dimension of the problem

        nodeDataAll = {} # all nodes in the gmsh model
        nodeTags, coord, parametricCoord = gmshModel.mesh.getNodes()
        numNodesAll = len(nodeTags)
        print("load node data")
        with alive_bar(numNodesAll) as bar:
            for ii in range(numNodesAll):
                nodeDataAll[nodeTags[ii]] = coord[3*ii:3*(ii+1)]
                bar()
        elementDataAll = {} # all element in the model
        elementTypes, elementTags, nodeTags = gmshModel.mesh.getElements(self.dimension)# element stores by its type
        nodeMap = {} #all active nodes, use dict to avoid duplicating
        eleMap = {}  #all active elements, use dict to avoid duplicating
        for ii in range(len(elementTypes)):
            eleType = elementTypes[ii]
            print("load element data: ",GmshModel.getType(eleType))
            elePros = gmshModel.mesh.getElementProperties(eleType)
            print(elePros)
            numberNodesPerElem = elePros[3]
            eleIds = elementTags[ii]
            nodeIds = nodeTags[ii]

            numEleType = len(eleIds)
            with alive_bar(numEleType) as bar:
                for jj in range(numEleType):
                    eid = eleIds[jj]
                    allNodes = nodeIds[numberNodesPerElem*jj:numberNodesPerElem*(jj+1)]
                    GmshModel.convertNodeNumbering(eleType,allNodes) # gmsh and oxfemm has not the same node ordering for some high-order elements
                    elementDataAll[eid] = allNodes
                    # update element map
                    if eid not in eleMap:
                        lastLength = len(eleMap)+1
                        eleMap[eid] = lastLength
                    # update node map
                    for ino in allNodes:
                        if ino not in nodeMap:
                            lastLength = len(nodeMap)+1
                            nodeMap[ino] = lastLength
                    bar()
        numNodes = len(nodeMap)
        self.nodes= np.zeros((numNodes,3))
        print(" node re-numbering")
        with alive_bar(numNodes) as bar:
            for key, val in nodeMap.items():
                self.nodes[val-1] =  nodeDataAll[key]
                bar()

        self.elementType = GmshModel.getType(elementTypes[0]) # get element type
        numElements= len(eleMap)
        self.elements = np.zeros((numElements,numberNodesPerElem),dtype=int)
        print(" node re-numbering in elements")
        with alive_bar(numEleType) as bar:
            for key, val in eleMap.items():
                for ii in range(len(elementDataAll[key])):
                    self.elements[val-1,ii] =  nodeMap[elementDataAll[key][ii]]
                bar()

        self.nodeMap = nodeMap
        self.elementMap = eleMap
        #get all groups, create node and element groups with a same name
        self.groupOfNodes = {}
        self.groupOfElements = {}
        self.groupOfFacets = {}
        #
        allPhysicalGroups = gmshModel.getPhysicalGroups()
        withFacets = False
        for phyGr in allPhysicalGroups:
            dimension = phyGr[0]
            phyIden  = phyGr[1]
            if dimension == self.dimension -1:
                withFacets = True
                break
        
        self.edges = np.zeros(0)
        if createEdges and self.dimension >1:
            allEdges = set()
            print("start building edges")
            with alive_bar(numElements) as bar:
                for key, values in elementDataAll.items():
                    # get edges in each element
                    allEdgesEls = GmshModel.getEdges(self.elementType,values)
                    numEdges = len(allEdgesEls)
                    # for each edge
                    for i in range(numEdges):
                        sube = [nodeMap[jj] for jj in allEdgesEls[i]]
                        sube.sort()
                        a = tuple(sube)
                        allEdges.add(a)
                    bar()
            print("done building edges, number of edges = ",len(allEdges))
            if len(allEdges)>0:
                firtsEdge = next(iter(allEdges))
                numberNodesPerEdge = len(firtsEdge)
                self.edges = np.zeros((len(allEdges),numberNodesPerEdge),dtype=int)
                with alive_bar(len(allEdges)) as bar:
                    lineId = 0
                    for ed in allEdges:
                        for ii in range(len(ed)):
                            self.edges[lineId,ii] = ed[ii]
                        lineId = lineId+1
                        bar()
            else:
                raise NotImplementedError("Three is no element in the current model")
                
        self.partitions = {}
        if nparts >1:
            print("start building partitions")
            self.gmshModel.mesh.partition(nparts)
            entities = self.gmshModel.getEntities()
            with alive_bar(numElements) as bar:
                for e in entities:
                    if e[0] == self.dimension:
                        procs = self.gmshModel.getPartitions(e[0], e[1])
                        elementTypes, elementTags, nodeTags =  self.gmshModel.mesh.getElements(e[0], e[1])
                        for par in procs:
                            if par not in self.partitions.keys():
                                self.partitions[par] = (set(),set())
                            nodeSet = self.partitions[par][0]
                            for nno in nodeTags:
                                for nn in nno:
                                    nodeSet.add(nodeMap[nn])
                            elSet = self.partitions[par][1]
                            for eleo in elementTags:
                                for ee in eleo:
                                    elSet.add(eleMap[ee])
                                    bar()
                
            print("done building partitions")
            total =0
            for par, val in self.partitions.items():
                print(f"part = {par} number of nodes = {len(val[0])} number of element = {len(val[1])}")
                total += len(val[1])
            print(f"total number of elements = {total}")
           
        self.internalFacets = {}
        self.externalFacets = {}
        if (createFacets or withFacets) and (self.dimension > 1):
            print("start building facets")
            with alive_bar(numElements) as bar:
                for key, values in elementDataAll.items():
                    # get facets in each element
                    allFacetEls = GmshModel.getFacets(self.elementType,values)
                    numFacets = len(allFacetEls)
                    # for each facet
                    for i in range(numFacets):
                        sube = [nodeMap[jj] for jj in allFacetEls[i]]
                        sube.sort()
                        a = tuple(sube)
                        # if facet exist already
                        if a in self.externalFacets.keys():
                            # get other
                            elem = self.externalFacets[a]
                            # remove from the list
                            self.externalFacets.pop(a)
                            # facet is an internal one
                            self.internalFacets[a] = (elem[0], (eleMap[key], i+1))
                        else:
                            self.externalFacets[a] = ((eleMap[key], i+1),)
                    bar()
            print("done building facets, number of facets = ",len(self.externalFacets)+len(self.internalFacets))
        
        
        for phyGr in allPhysicalGroups:
            dimension = phyGr[0]
            phyIden  = phyGr[1]
            # get nodes and element in this groups
            eleIds = set()
            nodeIds = set()
            allEntities = gmshModel.getEntitiesForPhysicalGroup(dimension, phyIden)
            for ent in allEntities:
                elementTypes, elementTags, nodeTags = gmshModel.mesh.getElements(dimension,ent)
                for ity in range(len(elementTypes)):
                    for ie in elementTags[ity]:
                        eleIds.add(ie)
                    for inod in nodeTags[ity]:
                        nodeIds.add(inod)

            # create a group of node
            self.groupOfNodes[(dimension,phyIden)] = np.array([nodeMap[ii] for ii in nodeIds])
            if dimension == self.dimension:
                # create a group of element
                self.groupOfElements[(dimension,phyIden)] = np.array([eleMap[ii] for ii in eleIds])

            elif dimension == self.dimension-1:
                # create a group of facets
                allFaceElements = {}
                for iele in eleIds:
                    # fine element
                    subelType, subnodeTags, subeldim, subeltag = gmshModel.mesh.getElement(iele)
                    subEleNodes = [nodeMap[k] for k in subnodeTags]
                    subEleNodes.sort()
                    a = tuple(subEleNodes)
                    if a in self.externalFacets.keys():
                        vals = self.externalFacets[a]
                        # take only one side when the internal facet is considered
                        eleIdex = vals[0][0]
                        position = vals[0][1]
                        if position in allFaceElements.keys():
                            allFaceElements[position].add(eleIdex)
                        else:
                            allFaceElements[position] = set()
                            allFaceElements[position].add(eleIdex)
                    else:
                        NotImplementedError(f"Element {iele} is not found in external facets")

                self.groupOfFacets[(dimension,phyIden)] = []
                for key, vals in allFaceElements.items():
                    self.groupOfFacets[(dimension,phyIden)].append([key, list(vals)])

            else:
                NotImplementedError(f"group with dimension {dimension} has not been implemented")

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
            allFacets={}
            internalFacets={}
            for ele in refineElements:
                allFacetEls = GmshModel.getFacets(self.elementType,self.elements[ele-1])
                numFacets = len(allFacetEls)
                # for each facet
                for i in range(numFacets):
                    sube = [jj for jj in allFacetEls[i]]
                    sube.sort()
                    a = tuple(sube)
                    # if facet exist already
                    if a in allFacets.keys():
                        # remove from the list
                        internalFacets[a] = (allFacets[a], [ele, i+1])
                        allFacets.pop(a)
                    else:
                        allFacets[a] = [ele, i+1]
            
            refineElementsOld = list(refineElements)
            refineElements.clear()
            
            #level1
            for key, val in allFacets.items():
                refineElements.append(val[0])
                if level >=1:
                    allFacetEls = GmshModel.getFacets(self.elementType,self.elements[val[0]-1])
                    numFacets = len(allFacetEls)
                    # for each facet
                    for i in range(numFacets):
                        if (i+1) != val[1]:
                            sube = [jj for jj in allFacetEls[i]]
                            sube.sort()
                            a = tuple(sube)
                            if a in internalFacets.keys():
                                ff=internalFacets[a]
                                refineElements.append(ff[0][0])
                                refineElements.append(ff[1][0])
                    
            #
            refineElements=[i for i in set(refineElements)]
        
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
        
        for key in self.groupOfFacets.keys():
            allPairs = []
            for ll in self.groupOfFacets[key]:
                for ele in ll[1]:
                    allPairs.append([ele,ll[0]])
            
            newPairs = {}
            for pp in allPairs:
                allFacets = self.getFacets(self.elementType, self.elements[pp[0]-1])
                facet = allFacets[pp[1]-1]
                facet.sort()
                facet = tuple(facet)
                found = False
                for newEle in newElementsMap[pp[0]]:
                    allFacetsNew = self.getFacets(self.elementType, newElements[newEle-1])
                    for ff in allFacetsNew:
                        ffT = list(ff)
                        ffT.sort()
                        ffT = tuple(ffT)
                        if facet == ffT:
                            if (allFacetsNew.index(ff)+1) not in newPairs.keys():
                                 newPairs[allFacetsNew.index(ff)+1] = []
                            newPairs[allFacetsNew.index(ff)+1].append(newEle)
                            found = True
                            break
                    if found:
                        break
                if not(found):
                    raise RuntimeError("something woring in facets reparing", pp)
            
            self.groupOfFacets[key] = [[k, v] for (k, v) in newPairs.items()]
            
        self.nodes = np.array(allNodes)
        self.elements = np.array(newElements,dtype=int)
     
        print("done mesh refining")
        


    def createPolyMesh(self, destFolder, boundaryConfig):
        """convert to polyMesh files for openFOAM
        Args:
            destFolder: location to save files
        Returns:
            None
        """
        curDir = os.getcwd()
        path = os.path.join(curDir,destFolder)
        if os.path.isdir(path):
            shutil.rmtree(path)
        os.makedirs(path)
        os.chdir(path)

        #write point
        ptFile = open("points","w")
        header = """/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2112                                  |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       vectorField;
    location    "constant/polyMesh";
    object      points;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n"""
        ptFile.write(header)
        ptFile.write("%d\n(\n"%len(self.nodes))
        for inod in range(len(self.nodes)):
            ptFile.write("(%.16g %.16g %.16g)\n"%(self.nodes[inod][0],self.nodes[inod][1],self.nodes[inod][2]))
        ptFile.write(")")
        ptFile.close()

        allBCs = []
        allBClist = set()
        for key, vals in self.groupOfFacets.items():
            data = set()
            for vv in vals:
                faceposition = vv[0]
                eleIds = vv[1]
                for ie in eleIds:
                    allSubs =GmshModel.getFacets(self.elementType,self.elements[ie-1])
                    subel = allSubs[faceposition-1]
                    subel.sort()
                    data.add(tuple(subel))
                    allBClist.add(tuple(subel))
            allBCs.append([key,data,0])
        allExternalBC = set(self.externalFacets.keys())
        allExternalBC.difference_update(allBClist)
        print("notDefined=",len(allExternalBC))

        # write faces
        ptFile = open("faces","w")
        header = """/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2112                                  |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       faceList;
    location    "constant/polyMesh";
    object      faces;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n"""
        ptFile.write(header)
        numFacets = len(self.externalFacets)+len(self.internalFacets)
        ptFile.write("%d\n(\n"%numFacets)
        owner = np.zeros(numFacets,dtype=int)
        neighbour=np.zeros(len(self.internalFacets),dtype=int)
        idex =0
        for face, elem in self.internalFacets.items():
            el1 = elem[0]
            el2 = elem[1]
            if el1[0] < el2[0]:
                owner[idex]=el1[0]-1
                pos = el1[1]-1
                neighbour[idex]=el2[0]-1
            else:
                owner[idex]=el2[0]-1
                pos=el2[1]-1
                neighbour[idex]=el1[0]-1
            allSubs =GmshModel.getFacets(self.elementType,self.elements[owner[idex]])
            numNod = len(allSubs[pos])
            ptFile.write("%d("%numNod)
            ptFile.write("%d"%(allSubs[pos][0]-1))
            for vv in range(1,numNod):
                ptFile.write(" %d"%(allSubs[pos][vv]-1))
            ptFile.write(")\n")
            idex = idex+1

        for vals in allBCs:
            vals[2] = idex
            for jj in vals[1]:
                elem = self.externalFacets[jj]
                owner[idex] =elem[0][0]-1
                pos=elem[0][1]-1
                allSubs =GmshModel.getFacets(self.elementType,self.elements[owner[idex]])
                numNod = len(allSubs[pos])
                ptFile.write("%d("%numNod)
                ptFile.write("%d"%(allSubs[pos][0]-1))
                for vv in range(1,numNod):
                    ptFile.write(" %d"%(allSubs[pos][vv]-1))
                ptFile.write(")\n")
                idex = idex+1
        for jj in allExternalBC:
            elem = self.externalFacets[jj]
            owner[idex] =elem[0][0]-1
            pos=elem[0][1]-1
            allSubs =GmshModel.getFacets(self.elementType,self.elements[owner[idex]])
            numNod = len(allSubs[pos])
            ptFile.write("%d("%numNod)
            ptFile.write("%d"%(allSubs[pos][0]-1))
            for vv in range(1,numNod):
                ptFile.write(" %d"%(allSubs[pos][vv]-1))
            ptFile.write(")\n")
            idex = idex+1

        ptFile.write(")")
        ptFile.close()

         #write owner
        ptFile = open("owner","w")
        header = """/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2112                                  |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    note        "nPoints:%d  nCells:%d  nFaces:%d  nInternalFaces:%d";
    class       labelList;
    location    "constant/polyMesh";
    object      owner;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n"""%(len(self.nodes),len(self.elements),numFacets,len(self.internalFacets))
        ptFile.write(header)
        ptFile.write("%d\n(\n"%len(owner))
        idex =0
        for ine in owner:
            ptFile.write("%d\n"%ine)
        ptFile.write(")")
        ptFile.close()

        #write neighbour
        ptFile = open("neighbour","w")
        header = """/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2112                                  |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    note        "nPoints:%d  nCells:%d  nFaces:%d  nInternalFaces:%d";
    class       labelList;
    location    "constant/polyMesh";
    object      neighbour;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n"""%(len(self.nodes),len(self.elements),numFacets,len(self.internalFacets))
        ptFile.write(header)
        ptFile.write("%d\n(\n"%len(neighbour))
        idex =0
        for ine in neighbour:
            ptFile.write("%d\n"%ine)
        ptFile.write(")")
        ptFile.close()

        #write neighbour
        ptFile = open("boundary","w")
        header = """/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2112                                  |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       polyBoundaryMesh;
    location    "constant/polyMesh";
    object      boundary;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n"""
        ptFile.write(header)
        ptFile.write("%d\n(\n"%len(allBCs))
        for ibc in allBCs:
            key = ibc[0]
            name = self.gmshModel.getPhysicalName(key[0],key[1])
            if name == "":
                name="BC"+str(key[1])
            if key[1] in boundaryConfig.keys():
                bcType = boundaryConfig[key[1]]
            elif name in boundaryConfig.keys():
                bcType = boundaryConfig[name]
            else:
                bcType = "patch"
            ptFile.write("\t%s\n\t{\n"%(name))
            if bcType=="patch":
                ptFile.write("\t\ttype\t\t%s;\n\t\tnFaces\t\t%d;\n\t\tstartFace\t\t%d;\n"%(bcType,len(ibc[1]),ibc[2]))
            else:
                ptFile.write("\t\ttype\t\t%s;\n\t\tinGroups\t\t1(%s);\n\t\tnFaces\t\t%d;\n\t\tstartFace\t\t%d;\n"%(bcType,bcType,len(ibc[1]),ibc[2]))
            ptFile.write("\t}\n\n")
        ptFile.write(")")
        ptFile.close()
        os.chdir(curDir)

