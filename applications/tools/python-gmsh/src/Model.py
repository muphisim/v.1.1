#class to define the mesh model
#Author(s): Van Dung NGUYEN, 2022
#

class Model:
    def getPartitions(self):
        """return all partitions
        Args:
            None
        Returns:
            {1: (nodeSet, elementSet), ...}
        """
        raise NotImplementedError('The method getPartitions() must be defined in Model!')
        
    def nodeInModel(self, nodeId):
        """return a boolean
        Args:
            nodeId: Node Id
        Returns:
            True if node is in the model, false otherwise
        """
        raise NotImplementedError('The method nodeInModel() must be defined in Model!')

    def elementInModel(self, elementId):
        """return a boolean
        Args:
            elementId: element Id
        Returns:
            True if element is in the model, false otherwise
        """
        raise NotImplementedError('The method elementInModel() must be defined in Model!')

    def getDimension(self):
        """return a dimension of the problem
        Args:
            None
            Returns:
                The dimension of the problem (2 or 3)
        """
        raise NotImplementedError('The method getDimension() must be defined in Model!')
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
        raise NotImplementedError('The method getAllNodes() must be defined in Model!')

    def getElementType(self):
        """return the element type of the mesh, only one element type is present
        Args:
            None
        Returns:
            The element type by string
        """
        raise NotImplementedError('The method getElementType() must be defined in Model!')

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
        raise NotImplementedError('The method getAllElements() must be defined in Model!')
        
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
        raise NotImplementedError('The method getAllEdges() must be defined in Model!')

    def getGroupOfElements(self, elementGroupIndex):
        """return a list of elements in a group identified elementGroupIndex
        Args:
            elementGroupIndex: unique number to identify the element group
        Returns:
            List of element identifiers.
        """
        raise NotImplementedError('The method getGroupOfElements() must be defined in Model!')

    def getGroupOfNodes(self, nodeGroupIndex):
        """return a list of nodes in a group identified nodeGroupIndex
        Args:
            nodeGroupIndex: unique number to identify the group of nodes
        Returns:
            List of node identifiers.
        """
        raise NotImplementedError('The method getGroupOfNodes() must be defined in Model!')

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
        raise NotImplementedError('The method getGroupOfFacets() must be defined in Model!')
        
    def createPolyMesh(self, destFolder):
        """convert to polyMesh files for openFOAM
        Args:
            destFolder: location to save files
        Returns:
            None
        """
        raise NotImplementedError('The method createPolyMesh() must be defined in Model!')
        
    def centroid_refine(self, boundaryOnly=False, elementGroups=None, level=1):
        """refine mesh form centroid
        Args:
            elementGroups: None if refinement is applied on all element; if not, the elements in this groups is refined
            boundaryOnly: True if only element close to the boundary is refined
        Returns:
            None
        """
        raise NotImplementedError('The method centroid_refine() must be defined in Model!')
