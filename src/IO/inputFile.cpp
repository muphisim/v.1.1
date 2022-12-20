//
//
// File authors: see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file inputFile.cpp
  \brief This file contains all functions related to the management of the input file (.inp). Please, be aware that most of the algorithms used here are not straightforward. If you are modifying, adding or deleting a function in this file (e.g., another constitutive model), it is strongly recommended to smartly copy/paste the code that already does the corresponding work.
*/
#include "inputFile.h"
#include "constitutiveList.h"
/*! \brief It reads the nodes in the domain
  @param[out] nodes Array of nodes
  @param[in] ndim Dimension of the problem
  @param[in] ifs Pointer to the input file (.inp)
  @param[in] npartShare Processor ownership of the nodes in the whole model npartShare[node i]=1 if the nodes is mine or shared with other processor, 0 otherwise (empty for sequential simulations)
*/
extern void nodesManagement(vector<classNodes *> &nodes, int &ndim, ifstream &ifs, vector<int> npartShare) {

    int maxsize = 100;
    classNodes *tempNode;
    bool flagNode = false;
    char tmpme[maxsize], tmpx[maxsize], tmpy[maxsize], tmpz[maxsize];
    memset(tmpme, 0, maxsize);
    memset(tmpx, 0, maxsize);
    memset(tmpy, 0, maxsize);
    memset(tmpz, 0, maxsize);
    string line;
    int tmp = 0, x = 0, y = 0, z = 0, ndimtmp = 0;
    int countNodLoc = 0;

    while (getline(ifs, line)) {
        if (line.compare(0, 5, "*Node") == 0) {
            flagNode = true;
        } else {
            if (line.compare(0, 1, "*") == 0 && flagNode) {
                break;
            } else {
                if (flagNode) {
                    tmp = int(line.size());
                    for (int j = 0; j < tmp; j++) {
                        if (line.compare(j, 1, ",") == 0) {
                            ndimtmp++;
                            if (x == 0) {
                                x = j;
                            } else if (y == 0) {
                                y = j;
                            } else {
                                z = j;
                            }
                        }
                    }
                    ndim = ndimtmp;
                    if (ndimtmp == 2) {
                        z = tmp + 1;
                    }
                    line.copy(tmpme, x, 0);
                    vector<double> xyz(3, 0);
                    if (npartShare.size() == 0) { // SEQUENTIAL
                        line.copy(tmpx, y - x - 1, x + 1);
                        line.copy(tmpy, z - y - 1, y + 1);
                        if (ndimtmp == 3) {
                            line.copy(tmpz, tmp - z, z + 1);
                            xyz[0] = atof(tmpx);
                            xyz[1] = atof(tmpy);
                            xyz[2] = atof(tmpz);
                        } else {
                            xyz[0] = atof(tmpx);
                            xyz[1] = atof(tmpy);
                            xyz[2] = 0;
                        }
                        tempNode = new classNodes(ndim, atoi(tmpme) - 1, xyz);
                        nodes.push_back(tempNode);
                    } else { // PARALLEL
                        if (npartShare[atoi(tmpme) - 1] == 1) {
                            line.copy(tmpx, y - x - 1, x + 1);
                            line.copy(tmpy, z - y - 1, y + 1);
                            if (ndimtmp == 3) {
                                line.copy(tmpz, tmp - z, z + 1);
                                xyz[0] = atof(tmpx);
                                xyz[1] = atof(tmpy);
                                xyz[2] = atof(tmpz);
                            } else {
                                xyz[0] = atof(tmpx);
                                xyz[1] = atof(tmpy);
                                xyz[2] = 0;
                            }
//                            tempNode= new classNodes(ndim, atoi(tmpme)-1, xyz);
                            tempNode = new classNodes(ndim, countNodLoc, xyz);
                            nodes.push_back(tempNode);
                            countNodLoc++;
                        }
                    }
                    memset(tmpme, 0, maxsize);
                    memset(tmpx, 0, maxsize);
                    memset(tmpy, 0, maxsize);
                    memset(tmpz, 0, maxsize);
                    x = 0;
                    y = 0;
                    z = 0;
                    ndimtmp = 0;
                }
            }
        }
    }
}

/*! \brief It reads all elements from the input file
  @param[out] elements Array with all elements in the domain
  @param[in] ndim Dimension of the problem
  @param[in] ifs Pointer to the input file (.inp)
  @param[in] nodes Array with all nodes in the domain
  @param[out] epart Processor ownership mapping of the elements in the whole model epart[element i]=rank (empty for sequential simulations)
  @param[in] nodLocal Mapping of local nodes to global nodes numbering nodlocal[local node index]= global node index (empty for sequential simulations)
*/
extern void elementsManagement(vector<classElements *> &elements, int ndim, ifstream &ifs, vector<classNodes *> &nodes,
                               vector<int> &epart, vector<int> &nodLocal, int &surface_order) {

    int maxsize = 100;
    bool flagElement = false;
    char tmpme2[maxsize], tempNode[maxsize];
    memset(tempNode, 0, maxsize);
    memset(tmpme2, 0, maxsize);

    int rank = 0;
#ifdef PARALLEL
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    CHKERRV(ierr);
#endif
    string line, type;
    int tmp = 0, tmp2 = 0, numNodes = 0;
    vector<int> myNodes;
    vector<int> commas;
    vector<int>::iterator it;
    int index = 0;
    classElements *tempElement;
    while (getline(ifs, line)) {
        if (line.compare(0, 8, "*Element") == 0) {
            tmp2 = int(line.size());
            if (line.compare(15, tmp2, "C3D4") == 0) {
                type.assign("C3D4");
                numNodes = 4; // For such type of element
                if (rank == 0) {
                    INFO("Type of elements %s", type.c_str());
                }
            } else if (line.compare(15, tmp2, "CPS3") == 0) {
                type.assign("CPS3");
                numNodes = 3; // For such type of element
                if (rank == 0) {
                    INFO("Type of elements %s", type.c_str());
                }
            } else if (line.compare(15, tmp2, "CPE4") == 0) {
                type.assign("CPE4");
                numNodes = 4; // For such type of element
                if (rank == 0) {
                    INFO("Type of elements %s", type.c_str());
                }
            } else if (line.compare(15, tmp2, "C3D10") == 0) {
                surface_order = 2;
                type.assign("C3D10");
                numNodes = 10; // For such type of element
                if (rank == 0) {
                    INFO("Type of elements %s", type.c_str());
                }
            } else if (line.compare(15, tmp2, "C3D8") == 0) {
                surface_order = 2;
                type.assign("C3D8");
                numNodes = 8; // For such type of element
                if (rank == 0) {
                    INFO("Type of elements %s", type.c_str());
                }
            } else if (line.compare(15, tmp2, "C3D10_lin_surf") == 0) {
                surface_order = 1;
                type.assign("C3D10");
                numNodes = 10; // For such type of element
                if (rank == 0) {
                    INFO("Type of elements %s", type.c_str());
                }
            } else if (line.compare(15, tmp2, "CPS6") == 0) {
                surface_order = 2;
                type.assign("CPS6");
                numNodes = 6; // For such type of element
                if (rank == 0) {
                    INFO("Type of elements %s", type.c_str());
                }
            } else if (line.compare(15, tmp2, "CPS6_lin_surf") == 0) {
                surface_order = 1;
                type.assign("CPS6");
                numNodes = 6; // For such type of element
                if (rank == 0) {
                    INFO("Type of elements %s", type.c_str());
                }
            } else {
                if (rank == 0) {
                    ERROR("Such type of element is not implemented in MuPhiSim %s", line.c_str());;
                    exit(-1);
                }
            }
            flagElement = true;
        } else {
            if (line.compare(0, 1, "*") == 0 && flagElement) {
                break;
            } else {
                if (flagElement) {
                    tmp = int(line.size());
                    int contc = 0;
                    for (int j = 0; j < tmp; j++) {
                        if (line.compare(j, 1, ",") == 0) {
                            commas.push_back(j);
                            contc++;
                        }
                    }
                    line.copy(tmpme2, commas[0], 0);
                    if (nodLocal.size() == 0) { // Sequential
                        for (int k = 1; k < numNodes; k++) {
                            memset(tempNode, 0, maxsize);
                            line.copy(tempNode, commas[k] - commas[k - 1] - 1, commas[k - 1] + 1);
                            myNodes.push_back(atoi(tempNode) - 1);
                            memset(tempNode, 0, maxsize);
                        }
                        line.copy(tempNode, tmp - commas[numNodes - 1] - 1, commas[numNodes - 1] + 1);
                        myNodes.push_back(atoi(tempNode) - 1);
                        if (type == "C3D4") {
                            tempElement = new classLinearTetrahedra(atoi(tmpme2) - 1, type, myNodes, ndim, nodes,
                                                                    atoi(tmpme2) - 1);
                        } else if (type == "C3D10") {
                            tempElement = new classQuadTetrahedra(atoi(tmpme2) - 1, type, myNodes, ndim, nodes,
                                                                  atoi(tmpme2) - 1);
                        } else if (type == "C3D8") {
                            tempElement = new classLinearHexahedra(atoi(tmpme2) - 1, type, myNodes, ndim, nodes,
                                                                  atoi(tmpme2) - 1);
                        } else if (type == "C3D10_lin_surf") {
                            tempElement = new classQuadTetrahedra(atoi(tmpme2) - 1, type, myNodes, ndim, nodes,
                                                                  atoi(tmpme2) - 1);
                        } else if (type == "CPS3") {
                            tempElement = new classLinearTriangles(atoi(tmpme2) - 1, type, myNodes, ndim, nodes,
                                                                   atoi(tmpme2) - 1);
                        } else if (type == "CPE4") {
                            tempElement = new classLinearSquare(atoi(tmpme2) - 1, type, myNodes, ndim, nodes,
                                                                   atoi(tmpme2) - 1);
                        } else if (type == "CPS6") {
                            tempElement = new classQuadTriangles(atoi(tmpme2) - 1, type, myNodes, ndim, nodes,
                                                                 atoi(tmpme2) - 1);
                        } else if (type == "CPS6_lin_surf") {
                            tempElement = new classQuadTriangles(atoi(tmpme2) - 1, type, myNodes, ndim, nodes,
                                                                 atoi(tmpme2) - 1);
                        } else if (type == "Segments") {
                            tempElement = new classLinearLine(atoi(tmpme2) - 1, type, myNodes, ndim, nodes,
                                                              atoi(tmpme2) - 1);
                        } else {
                            ERROR("The type of element %s is not available in MuPhiSim", type.c_str());
                            exit(-1);
                        }
                        elements.push_back(tempElement);
                    } else { //parallel
                        if (epart[atoi(tmpme2) - 1] == rank) {
                            for (int k = 1; k < numNodes; k++) {
                                memset(tempNode, 0, maxsize);
                                line.copy(tempNode, commas[k] - commas[k - 1] - 1, commas[k - 1] + 1);
                                myNodes.push_back(atoi(tempNode) - 1);//XXX origins for neighbours
                                memset(tempNode, 0, maxsize);
                            }
                            line.copy(tempNode, tmp - commas[numNodes - 1] - 1, commas[numNodes - 1] + 1);
                            myNodes.push_back(atoi(tempNode) - 1);


                            for (int k = 0; k < numNodes; k++) {
                                it = find(nodLocal.begin(), nodLocal.end(), myNodes[k]);
                                myNodes[k] = it - nodLocal.begin();
                            }

                            if (type == "C3D4") {
                                tempElement = new classLinearTetrahedra(elements.size(), type, myNodes, ndim, nodes,
                                                                        elements.size());
                            } else if (type == "C3D10") {
                                tempElement = new classQuadTetrahedra(elements.size(), type, myNodes, ndim, nodes,
                                                                      elements.size());
                            } else if (type == "C3D8") {
                                tempElement = new classLinearHexahedra(elements.size(), type, myNodes, ndim, nodes,
                                                                      elements.size());
                            } else if (type == "C3D10_lin_surf") {
                                tempElement = new classQuadTetrahedra(elements.size(), type, myNodes, ndim, nodes,
                                                                      elements.size());
                            } else if (type == "CPS3") {

                                tempElement = new classLinearTriangles(elements.size(), type, myNodes, ndim, nodes,
                                                                       elements.size());
                            } else if (type == "CPE4") {

                                tempElement = new classLinearSquare(elements.size(), type, myNodes, ndim, nodes,
                                                                       elements.size());
                            } else if (type == "CPS6") {

                                tempElement = new classQuadTriangles(elements.size(), type, myNodes, ndim, nodes,
                                                                     elements.size());
                            } else if (type == "CPS6_lin_surf") {

                                tempElement = new classQuadTriangles(elements.size(), type, myNodes, ndim, nodes,
                                                                     elements.size());
                            } else {
                                ERROR("The type of element %s is not available in MuPhiSim", type.c_str());
                                exit(-1);
                            }
                            elements.push_back(tempElement);
                            index++;
                        }
                    }
                    memset(tmpme2, 0, maxsize);
                    commas.clear();
                    myNodes.clear();
                }
            }
        }
    }

    if (!flagElement) {
        ERROR("There is not a list of elements in the input file");
        exit(-1);
    }
}

/*! \brief It defines what nodes in the domain correspond to FEM nodes
  @param[out] nodesFEM Array of integers with all labels of the general array of nodes that are FEM nodes
  @param[in] ifs Pointer to the input file (.inp)
  @param[out]  nodes Array with all nodes in the domain
  @param[out] npartShare Processor ownership of the nodes in the whole model npartShare[node i]=1 if the nodes is mine or shared with other processor, 0 otherwise (empty for sequential simulations)
*/
extern void
nodesFEMManagement(vector<int> &nodesFEM, ifstream &ifs, vector<classNodes *> &nodes, vector<int> &npartShare) {

    int maxsize = 100;
    char tempNode[maxsize];
    memset(tempNode, 0, maxsize);

    bool flagNode = false;
    string line;
    int tmp = 0, size = 0;

    while (getline(ifs, line)) {
        // cout<<line<<endl;
        if (line.compare(0, 10, "*FEM Nodes") == 0) {
            flagNode = true;
        } else {
            if (line.compare(0, 1, "*") == 0 && flagNode) {
                break;
            } else {
                if (line.compare(0, 3, "All") == 0 && flagNode) {
                    int nNod = nodes.size();
                    for (int i = 0; i < nNod; i++) {
                        nodesFEM.push_back(nodes[i]->getMe());
                    }
                    break;
                } else {
                    if (flagNode) {
                        size = int(line.size());
                        line.copy(tempNode, size, 0);
                        tmp = atoi(tempNode) - 1;
                        if (npartShare.size() == 0) {
                            nodesFEM.push_back(tmp);
                        } else {
                            if (npartShare[tmp] == 1) {
                                nodesFEM.push_back(tmp);
                            }
                        }
                        memset(tempNode, 0, maxsize);
                    }
                }
            }
        }
    }
    if (!flagNode) {
        WARNING("There is not a list of FEM nodes in the input file");
    }
}

extern void
nodesFEMManagement(vector<int> &nodesFEM, ifstream &ifs, vector<classNodes *> &nodes, vector<int> &npartShare,
                   vector<int> nodLocal) {

    int maxsize = 100;
    char tempNode[maxsize];
    memset(tempNode, 0, maxsize);

    bool flagNode = false;
    string line;
    int tmp = 0, size = 0;

    while (getline(ifs, line)) {
        if (line.compare(0, 10, "*FEM Nodes") == 0) {
            flagNode = true;
        } else {
            if (line.compare(0, 1, "*") == 0 && flagNode) {
                break;
            } else {
                if (line.compare(0, 3, "All") == 0 && flagNode) {
                    int nNod = nodes.size();
                    for (int i = 0; i < nNod; i++) {
                        nodesFEM.push_back(nodes[i]->getMe());
                    }
                    break;
                } else {
                    if (flagNode) {
                        size = int(line.size());
                        line.copy(tempNode, size, 0);
                        tmp = atoi(tempNode) - 1;
                        if (npartShare.size() == 0) {
                            nodesFEM.push_back(tmp);
                        } else {
                            if (npartShare[tmp] == 1) {
                                int nodeInLocal;
                                std::vector<int>::iterator it = std::find(nodLocal.begin(), nodLocal.end(), tmp);
                                if (it != nodLocal.end()) {
                                    nodeInLocal = it - nodLocal.begin();
                                }
                                nodesFEM.push_back(nodeInLocal);
                            }
                        }
                        memset(tempNode, 0, maxsize);
                    }
                }
            }
        }
    }
    if (!flagNode) {
        WARNING("There is not a list of FEM nodes in the input file");
    }
}

/*! \brief It defines what nodes in the domain correspond to MM nodes
  @param[out] nodesMM Array of integers with all labels of the general array of nodes that are MM nodes
  @param[in] ifs Pointer to the input file (.inp)
  @param[in]  nodes Array with all nodes in the domain
  @param[in] npartShare Processor ownership of the nodes in the whole model npartShare[node i]=1 if the nodes is mine or shared with other processor, 0 otherwise (empty for sequential simulations)
*/
extern void
nodesMMManagement(vector<int> &nodesMM, ifstream &ifs, vector<classNodes *> &nodes, vector<int> &npartShare) {

    int maxsize = 100;
    char tempNode[maxsize];
    memset(tempNode, 0, maxsize);
    bool flagNode = false;
    string line;
    int tmp = 0, size = 0;

    while (getline(ifs, line)) {

        if (line.compare(0, 9, "*MM Nodes") == 0) {
            flagNode = true;
        } else {
            if (line.compare(0, 1, "*") == 0 && flagNode) {
                break;
            } else {
                if (line.compare(0, 3, "All") == 0 && flagNode) {
                    int nNod = nodes.size();
                    for (int i = 0; i < nNod; i++) {
                        nodesMM.push_back(nodes[i]->getMe());
                    }
                    break;
                } else {
                    if (flagNode) {
                        size = int(line.size());
                        line.copy(tempNode, size, 0);
                        tmp = atoi(tempNode) - 1;
                        if (npartShare.size() == 0) {
                            nodesMM.push_back(tmp);
                        } else {
                            if (npartShare[tmp] == 1) {
                                nodesMM.push_back(tmp);
                            }
                        }
                        memset(tempNode, 0, maxsize);
                    }
                }
            }
        }
    }
    if (!flagNode) {
        WARNING("There is not a list of MM nodes in the input file");
    }
}

extern void nodesMMManagement(vector<int> &nodesMM, ifstream &ifs, vector<classNodes *> &nodes, vector<int> &npartShare,
                              vector<int> nodLocal) {

    int maxsize = 100;
    char tempNode[maxsize];
    memset(tempNode, 0, maxsize);
    bool flagNode = false;
    string line;
    int tmp = 0, size = 0;

    while (getline(ifs, line)) {
        if (line.compare(0, 9, "*MM Nodes") == 0) {
            flagNode = true;
        } else {
            if (line.compare(0, 1, "*") == 0 && flagNode) {
                break;
            } else {
                if (line.compare(0, 3, "All") == 0 && flagNode) {
                    int nNod = nodes.size();
                    for (int i = 0; i < nNod; i++) {
                        nodesMM.push_back(nodes[i]->getMe());
                    }
                    break;
                } else {
                    if (flagNode) {
                        size = int(line.size());
                        line.copy(tempNode, size, 0);
                        tmp = atoi(tempNode) - 1;
                        if (npartShare.size() == 0) {
                            nodesMM.push_back(tmp);
                        } else {
                            if (npartShare[tmp] == 1) {
                                int nodeInLocal;

                                std::vector<int>::iterator it = std::find(nodLocal.begin(), nodLocal.end(), tmp);
                                if (it != nodLocal.end()) {
                                    nodeInLocal = it - nodLocal.begin();
                                }
                                nodesMM.push_back(nodeInLocal);
                            }
                        }
                        memset(tempNode, 0, maxsize);
                    }
                }
            }
        }
    }
    if (!flagNode) {

        WARNING("There is not a list of MM nodes in the input file");
    }
}


/*! \brief It reads the shared nodes between FEM and MM domains
  @param[out] nodesShare Array of integers with all labels of the general array of nodes that are shared nodes
  @param[in] ifs Pointer to the input file (.inp)
  @param[out] npartShare Processor ownership of the nodes in the whole model npartShare[node i]=1 if the nodes is mine or shared with other processor, 0 otherwise (empty for sequential simulations)
*/
extern void nodesShareManagement(vector<int> &nodesShare, ifstream &ifs, vector<int> &npartShare) {

    int maxsize = 100;
    char tempNode[maxsize];
    memset(tempNode, 0, maxsize);

    bool flagNode = false;
    string line;
    int tmp = 0, size = 0;
    while (getline(ifs, line)) {
        if (line.compare(0, 13, "*Shared Nodes") == 0) {
            flagNode = true;
        } else {
            if (line.compare(0, 1, "*") == 0 && flagNode) {
                break;
            } else {
                if (flagNode) {
                    size = int(line.size());
                    line.copy(tempNode, size, 0);
                    tmp = atoi(tempNode) - 1;
                    if (npartShare.size() == 0) {
                        nodesShare.push_back(tmp);
                    } else {
                        if (npartShare[tmp] == 1) {
                            nodesShare.push_back(tmp);
                        }
                    }
                    memset(tempNode, 0, maxsize);
                }
            }
        }
    }
    if (!flagNode) {

        WARNING("There is not a list of Shared nodes in the input file");
    }
}

/*! \brief It defines what elements in the domain that correspond to FEM elements
  @param[out] elementsFEM Array of integers with all labels of the general array of elements that are FEM elements
  @param[in] ifs Pointer to the input file (.inp)
  @param[in]  elements Array with all elements (or background cells) in the domain
  @param[in] epart Processor ownership mapping of the elements in the whole model epart[element i]=rank (empty for sequential simulations)
*/
extern void
elementsFEMManagement(vector<int> &elementsFEM, ifstream &ifs, vector<classElements *> &elements, vector<int> &epart) {

    int maxsize = 100;
    char tempNode[maxsize];
    memset(tempNode, 0, maxsize);
    bool flagNode = false;
    string line;
    int tmp = 0, size = 0;


    int rank = 0;
#ifdef PARALLEL
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    CHKERRV(ierr);
#endif

    int elemFEMCount = 0;
    while (getline(ifs, line)) {
        if (line.compare(0, 13, "*FEM Elements") == 0) {
            flagNode = true;
        } else {
            if (line.compare(0, 1, "*") == 0 && flagNode) {
                break;
            } else {
                if (line.compare(0, 3, "All") == 0 && flagNode) {
                    int nEle = elements.size();
                    for (int i = 0; i < nEle; i++) {
                        elementsFEM.push_back(elements[i]->getMe());
                    }
                    break;
                } else {
                    if (flagNode) {
                        size = int(line.size());
                        line.copy(tempNode, size, 0);
                        tmp = atoi(tempNode) - 1;
                        if (epart.size() == 0) {
                            elementsFEM.push_back(tmp);
                        } else {
                            if (epart[tmp] == rank) {
                                elementsFEM.push_back(elemFEMCount);
                                elemFEMCount++;
                            }
                        }
                        memset(tempNode, 0, maxsize);
                    }
                }
            }
        }
    }
    if (!flagNode) {
        WARNING("There is not a list of FEM elements in the input file");
    }

}

extern void elementsFEMManagement(vector<int> &elementsFEM, ifstream &ifs, vector<classElements *> &elements,
                                  vector<int> &epart, vector<int> elemLocal) {

    int maxsize = 100;
    char tempNode[maxsize];
    memset(tempNode, 0, maxsize);
    bool flagNode = false;
    string line;
    int tmp = 0, size = 0;


    int rank = 0;
#ifdef PARALLEL
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    CHKERRV(ierr);
#endif

    while (getline(ifs, line)) {
        if (line.compare(0, 13, "*FEM Elements") == 0) {
            flagNode = true;
        } else {
            if (line.compare(0, 1, "*") == 0 && flagNode) {
                break;
            } else {
                if (line.compare(0, 3, "All") == 0 && flagNode) {
                    int nEle = elements.size();
                    for (int i = 0; i < nEle; i++) {
                        elementsFEM.push_back(elements[i]->getMe());
                    }
                    break;
                } else {
                    if (flagNode) {
                        size = int(line.size());
                        line.copy(tempNode, size, 0);
                        tmp = atoi(tempNode) - 1;
                        if (epart.size() == 0) {
                            elementsFEM.push_back(tmp);
                        } else {
                            if (epart[tmp] == rank) {
                                int elemInLocal;

                                std::vector<int>::iterator it = std::find(elemLocal.begin(), elemLocal.end(), tmp);
                                if (it != elemLocal.end()) {
                                    elemInLocal = it - elemLocal.begin();
                                }
                                elementsFEM.push_back(elemInLocal);
                            }
                        }
                        memset(tempNode, 0, maxsize);
                    }
                }
            }
        }
    }
    if (!flagNode) {
        // if(sched_getcpu()==0)
        // {
        WARNING("There is not a list of FEM elements in the input file");
        // }
    }

}


/*! \brief It defines what elements in the domain that correspond to MM integration cells
  @param[out] elementsMM Array of integers with all labels of the general array of elements that are MM integration cells
  @param[in] ifs Pointer to the input file (.inp)
  @param[in]  elements Array with all elements (or background cells) in the domain
  @param[in] epart Processor ownership mapping of the elements in the whole model epart[element i]=rank (empty for sequential simulations)
*/
extern void elementsMMManagement(vector<int> &elementsMM, ifstream &ifs, vector<classElements *> &elements,
                                 vector<int> &epart) {

    int maxsize = 100;
    char tempNode[maxsize];
    memset(tempNode, 0, maxsize);

    bool flagNode = false;
    string line;
    int tmp = 0, size = 0;

    int rank = 0;
#ifdef PARALLEL
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    CHKERRV(ierr);
#endif
    int elemMMCount = 0;
    while (getline(ifs, line)) {
        if (line.compare(0, 12, "*MM Elements") == 0) {
            flagNode = true;
        } else {
            if (line.compare(0, 1, "*") == 0 && flagNode) {
                break;
            } else {
                if (line.compare(0, 3, "All") == 0 && flagNode) {
                    int nEle = elements.size();
                    for (int i = 0; i < nEle; i++) {
                        elementsMM.push_back(elements[i]->getMe());
                    }
                    break;
                } else {
                    if (flagNode) {
                        size = int(line.size());
                        line.copy(tempNode, size, 0);
                        tmp = atoi(tempNode) - 1;
                        if (epart.size() == 0) {
                            elementsMM.push_back(tmp);
                        } else {
                            if (epart[tmp] == rank) {
                                elementsMM.push_back(elemMMCount);
                                elemMMCount++;
                            }
                        }
                        memset(tempNode, 0, maxsize);
                    }
                }
            }
        }
    }
    if (!flagNode) {
        WARNING("There is not a list of MM elements in the input file");
    }
}

extern void elementsMMManagement(vector<int> &elementsMM, ifstream &ifs, vector<classElements *> &elements,
                                 vector<int> &epart, vector<int> elemLocal) {

    int maxsize = 100;
    char tempNode[maxsize];
    memset(tempNode, 0, maxsize);

    bool flagNode = false;
    string line;
    int tmp = 0, size = 0;

    int rank = 0;
#ifdef PARALLEL
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    CHKERRV(ierr);
#endif
    while (getline(ifs, line)) {
        if (line.compare(0, 12, "*MM Elements") == 0) {
            flagNode = true;
        } else {
            if (line.compare(0, 1, "*") == 0 && flagNode) {
                break;
            } else {
                if (line.compare(0, 3, "All") == 0 && flagNode) {
                    int nEle = elements.size();
                    for (int i = 0; i < nEle; i++) {
                        elementsMM.push_back(elements[i]->getMe());
                    }
                    break;
                } else {
                    if (flagNode) {
                        size = int(line.size());
                        line.copy(tempNode, size, 0);
                        tmp = atoi(tempNode) - 1;
                        if (epart.size() == 0) {
                            elementsMM.push_back(tmp);
                        } else {
                            if (epart[tmp] == rank) {
                                int elemInLocal;

                                std::vector<int>::iterator it = std::find(elemLocal.begin(), elemLocal.end(), tmp);
                                if (it != elemLocal.end()) {
                                    elemInLocal = it - elemLocal.begin();
                                }
                                elementsMM.push_back(elemInLocal);
                            }
                        }
                        memset(tempNode, 0, maxsize);
                    }
                }
            }
        }
    }
    if (!flagNode) {
        WARNING("There is not a list of MM elements in the input file");
    }
}

/*! \brief It assigns the constitutive model to the elements. Please note that different constitutive models can be assigned to different regions in the domain
  @param[inout] elements Array with all elements (or background cells) in the domain
  @param[in] ifs Pointer to the input file (.inp)
  @param[in] ndim Dimension of the problem
  @param[out] consMod Array with the constitutive models in the domain
  @param[in] epart Processor ownership mapping of the elements in the whole model epart[element i]=rank (empty for sequential simulations)
  @param[inout] elemLocal Mapping of local elements to global elements numbering -> elemlocal[local element index] = global element index (empty for sequential simulations)
*/
extern void materialManagement(vector<classElements *> &elements, ifstream &ifs, int ndim,
                               vector<constitutiveModels *> &consMod, vector<int> &epart, vector<int> &elemLocal) {

    int maxsize = 100;
    char tempNode[maxsize];
    memset(tempNode, 0, maxsize);
    bool flagNode = false, flagList = false;
    string line, templine;
    string name;
    string name_cons;
    int size = 0;
    double rho = 0;
    constitutiveModels *constitutiveModel = NULL;
    char tmpme2[maxsize];
    memset(tmpme2, 0, maxsize);
    int tmp = 0;
    vector<int> commas;
    vector<int>::iterator it;
    int rank = 0;
    bool flagExtra = false;
    bool flagStochastic = false;
    bool stochasticMapping = false;
    vector<int> listElements;
#ifdef PARALLEL
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    CHKERRV(ierr);
#endif

    while (getline(ifs, line)) {
        if (line.compare(0, 9, "*MATERIAL") == 0) {
            while (!consMod.size() == 0) { //consMod should be empty...
                consMod.pop_back(); //otherwise we need to delete the cons models of previous region before defining a new one
            }
            listElements.clear(); //the list of elements is cleaned for next region
            flagNode = true;
            flagList = false;
            size = int(line.size());
            name.assign(line, 16, size - 16);
            size = 0;
        } else {
            if (line.compare(0, 8, "*Density") == 0 && flagNode) {
                getline(ifs, templine);
                size = int(templine.size());
                templine.copy(tempNode, size, 0);
                rho = atof(tempNode);
                memset(tempNode, 0, maxsize);
                templine.clear();
            } else if (line.compare(0, 15, "*Linear-Elastic") == 0 && flagNode) {
                getline(ifs, templine);
                tmp = int(templine.size());
                int contc = 0;
                commas.clear();
                for (int j = 0; j < tmp; j++) {
                    if (templine.compare(j, 1, ",") == 0) {
                        commas.push_back(j);
                        contc++;
                    }
                }
                commas.push_back(tmp + 1);
                templine.copy(tmpme2, commas[0], 0);
                double young = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                templine.copy(tempNode, commas[1] - commas[0] - 1, commas[0] + 1);
                double nu = atof(tempNode);
                memset(tempNode, 0, maxsize);
                int stressState = 0;//default plane stress
                if (commas.size() == 3) 
                {
                    templine.copy(tempNode, commas[2] - commas[1] - 1, commas[1] + 1);
                    string stressStateString = "";
                    for (int i = 0; i < sizeof(tempNode) / sizeof(char); i++) {
                        if (tempNode[i] != ' ')//remove all whitespaces
                        {
                            stressStateString = stressStateString + tempNode[i];
                        }
                    }
                    if (stressStateString.compare(0, 11, "PlaneStress") == 0) {
                        stressState = 0;
                    } else if (stressStateString.compare(0, 11, "PlaneStrain") == 0) {
                        stressState = 1;
                    } else {
                        stressState = 0;
                    }
                    memset(tempNode, 0, maxsize);
                }
                name_cons = "Linear-Elastic";
                constitutiveModel = new linearElastic(ndim, name, name_cons, rho, young, nu,
                                                      stressState);
                commas.clear();
                consMod.push_back(constitutiveModel);
            } else if (line.compare(0, 13, "*HyperElastic") == 0 && flagNode) {
                if ((!line.compare(15, 26, "Neo-Hookean") == 0) &&
                    (!line.compare(15, 43, "3D-Compressible-Neo-Hookean") == 0) &&
                    (!line.compare(15, 34, "St-Venant-Kirchhoff") == 0) &&
                    (!line.compare(15, 49, "St-Venant-Kirchhoff-ThermoMechanics") == 0)) {
                    ERROR("The %s model is not implemented in MuPhiSim", line.c_str());
                    exit(-1);
                }
                getline(ifs, templine);
                tmp = int(templine.size());
                int contc = 0;
                for (int j = 0; j < tmp; j++) {
                    if (templine.compare(j, 1, ",") == 0) {
                        commas.push_back(j);
                        contc++;
                    }
                }
                commas.push_back(tmp + 1);
                templine.copy(tmpme2, commas[0], 0);
                double young = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                templine.copy(tempNode, commas[1] - commas[0] - 1, commas[0] + 1);
                double nu = atof(tempNode);
                memset(tempNode, 0, maxsize);
                int stressState = 0;//default plane stress
                double alpha;
                double theta_ini;
                if (line.compare(15, 49, "St-Venant-Kirchhoff-ThermoMechanics") == 0) {
                    templine.copy(tempNode, commas[2] - commas[1] - 1, commas[1] + 1);
                    alpha = atof(tempNode);
                    memset(tempNode, 0, maxsize);
                    templine.copy(tempNode, commas[3] - commas[2] - 1, commas[2] + 1);
                    theta_ini = atof(tempNode);
                    memset(tempNode, 0, maxsize);
                    name_cons = "St-Venant-Kirchhoff-ThermoMechanics";

                    if (commas.size() == 5)// non default
                    {
                        templine.copy(tempNode, commas[4] - commas[3] - 1, commas[3] + 1);
                        string stressStateString = "";

                        for (int i = 0; i < sizeof(tempNode) / sizeof(char); i++) {
                            if (tempNode[i] != ' ')//remove all whitespaces
                            {
                                stressStateString = stressStateString + tempNode[i];
                            }
                        }
                        if (stressStateString.compare(0, 11, "PlaneStress") == 0) {
                            stressState = 0;
                        } else if (stressStateString.compare(0, 11, "PlaneStrain") == 0) {
                            stressState = 1;
                        } else {
                            stressState = 0;
                        }
                        memset(tempNode, 0, maxsize);
                    }
                } else {
                    if (commas.size() == 3)// non default
                    {
                        templine.copy(tempNode, commas[2] - commas[1] - 1, commas[1] + 1);
                        string stressStateString = "";

                        for (int i = 0; i < sizeof(tempNode) / sizeof(char); i++) {
                            if (tempNode[i] != ' ')//remove all whitespaces
                            {
                                stressStateString = stressStateString + tempNode[i];
                            }
                        }
                        if (stressStateString.compare(0, 11, "PlaneStress") == 0) {
                            stressState = 0;
                        } else if (stressStateString.compare(0, 11, "PlaneStrain") == 0) {
                            stressState = 1;
                        } else {
                            stressState = 0;
                        }
                        memset(tempNode, 0, maxsize);
                    }
                }

                if (line.compare(15, 26, "Neo-Hookean") == 0) {
                    name_cons = "Neo-Hookean";
                    constitutiveModel = new classHyperElasticNeoHookean(ndim, name, name_cons, rho, young, nu);
                    commas.clear();
                } else if (line.compare(15, 43, "3D-Compressible-Neo-Hookean") == 0) {
                    name_cons = "3D-Compressible-Neo-Hookean";
                    constitutiveModel = new classHyperElastic3DCompressibleNeoHookean(ndim, name, name_cons, rho, young, nu);
                    commas.clear();
                } else if (line.compare(15, 34, "St-Venant-Kirchhoff") == 0) {
                    name_cons = "St-Venant-Kirchhoff";
                    constitutiveModel = new classHyperElasticStVenantKirchhoff(ndim, name, name_cons, rho, young, nu,
                                                                               stressState);
                    commas.clear();
                } else if (line.compare(15, 50, "St-Venant-Kirchhoff-ThermoMechanics") == 0) {
                    constitutiveModel = new classHyperElasticStVenantKirchhoffThermoMechanics(ndim, name, name_cons,
                                                                                              rho,
                                                                                              young, nu, stressState,
                                                                                              alpha, theta_ini);
                    commas.clear();

                }
                consMod.push_back(constitutiveModel);
            } else if (line.compare(0, 48, "*3DPrinting-St-Venant-Kirchhoff-Thermoelasticity") == 0 && flagNode) {
                templine.clear();
                getline(ifs, templine);
                tmp = int(templine.size());
                commas.clear();
                int contc = 0;
                for (int j = 0; j < tmp; j++) {
                    if (templine.compare(j, 1, ",") == 0) {
                        commas.push_back(j);
                        contc++;
                    }
                }
                commas.push_back(tmp + 1);
                templine.copy(tmpme2, commas[0], 0);
                double rhoLiquid = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                templine.copy(tmpme2, commas[1] - commas[0] - 1, commas[0] + 1);
                double rho_PL = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                templine.copy(tmpme2, commas[2] - commas[1] - 1, commas[1] + 1);
                double youngRefPowder = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                templine.copy(tmpme2, commas[3] - commas[2] - 1, commas[2] + 1);
                double young1_PL = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                templine.copy(tmpme2, commas[4] - commas[3] - 1, commas[3] + 1);
                double youngRefSolid = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                templine.copy(tmpme2, commas[5] - commas[4] - 1, commas[4] + 1);
                double young1_SL = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                templine.copy(tmpme2, commas[6] - commas[5] - 1, commas[5] + 1);
                double young2_SL = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                templine.copy(tmpme2, commas[7] - commas[6] - 1, commas[6] + 1);
                double nuRefSolid = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                templine.copy(tmpme2, commas[8] - commas[7] - 1, commas[7] + 1);
                double nu_SL = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                templine.copy(tmpme2, commas[9] - commas[8] - 1, commas[8] + 1);
                double  alphaRefSolid = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                templine.copy(tmpme2, commas[10] - commas[9] - 1, commas[9] + 1);
                double logGrowthRate = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                templine.copy(tmpme2, commas[11] - commas[10] - 1, commas[10] + 1);
                double  theta_sigMidpoint = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                templine.copy(tmpme2, commas[12] - commas[11] - 1, commas[11] + 1);
                double theta_ref = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                templine.copy(tmpme2, commas[13] - commas[12] - 1, commas[12] + 1);
                double theta_SL = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                templine.copy(tmpme2, commas[14] - commas[13] - 1, commas[13] + 1);
                double theta_liq = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                commas.clear();
                name_cons = "3DprintingThermoelasticity";
		constitutiveModel = new class3DprintingThermoelasticityStVenantKirchhoff(ndim, name, name_cons, rho,
                                                                                         rhoLiquid, rho_PL,
                                                                                         youngRefPowder, young1_PL,
                                                                                         youngRefSolid, young1_SL,
                                                                                         young2_SL, nuRefSolid,
                                                                                         nu_SL, alphaRefSolid,
                                                                                         logGrowthRate,
                                                                                         theta_sigMidpoint, 
                                                                                         theta_ref,
                                                                                         theta_SL, theta_liq);
                consMod.push_back(constitutiveModel);
            } else if (line.compare(0, 16, "*IsoMorphoGrowth") == 0 && flagNode) {
                templine.clear();
                getline(ifs, templine);
                tmp = int(templine.size());
                commas.clear();
                int contc = 0;
                for (int j = 0; j < tmp; j++) {
                    if (templine.compare(j, 1, ",") == 0) {
                        commas.push_back(j);
                        contc++;
                    }
                }
                commas.push_back(tmp + 1);
                templine.copy(tmpme2, commas[0], 0);
                double young = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                templine.copy(tempNode, commas[1] - commas[0] - 1, commas[0] + 1);
                double nu = atof(tempNode);
                memset(tempNode, 0, maxsize);
                templine.copy(tempNode, commas[2] - commas[1] - 1, commas[1] + 1);
                double Gc = atof(tempNode);
                memset(tempNode, 0, maxsize);
                commas.clear();
                name_cons = "IsoMorphoGrowth";
                constitutiveModel = new classIsoMorphogeneticGrowth(ndim, name, name_cons, rho, young, nu, Gc);
                consMod.push_back(constitutiveModel);
            } else if (line.compare(0, 12, "*FiberGrowth") == 0 && flagNode) {
                getline(ifs, templine);
                tmp = int(templine.size());
                int contc = 0;
                for (int j = 0; j < tmp; j++) {
                    if (templine.compare(j, 1, ",") == 0) {
                        commas.push_back(j);
                        contc++;
                    }
                }
                commas.push_back(tmp + 1);
                templine.copy(tmpme2, commas[0], 0);
                double young = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                templine.copy(tempNode, commas[1] - commas[0] - 1, commas[0] + 1);
                double nu = atof(tempNode);
                memset(tempNode, 0, maxsize);
                templine.copy(tempNode, commas[2] - commas[1] - 1, commas[1] + 1);
                double Gc = atof(tempNode);
                memset(tempNode, 0, maxsize);
                commas.clear();
                name_cons = "FiberGrowth";
                constitutiveModel = new classFiberGrowth(ndim, name, name_cons, rho, young, nu, Gc);
                consMod.push_back(constitutiveModel);
            } else if (line.compare(0, 11, "*AreaGrowth") == 0 && flagNode) {
                getline(ifs, templine);
                tmp = int(templine.size());
                int contc = 0;
                for (int j = 0; j < tmp; j++) {
                    if (templine.compare(j, 1, ",") == 0) {
                        commas.push_back(j);
                        contc++;
                    }
                }
                commas[2] = tmp + 1;
                templine.copy(tmpme2, commas[0], 0);
                double young = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                templine.copy(tempNode, commas[1] - commas[0] - 1, commas[0] + 1);
                double nu = atof(tempNode);
                memset(tempNode, 0, maxsize);
                templine.copy(tempNode, commas[2] - commas[1] - 1, commas[1] + 1);
                double Gc = atof(tempNode);
                memset(tempNode, 0, maxsize);
                commas.clear();
                name_cons = "AreaGrowth";
                constitutiveModel = new classAreaGrowth(ndim, name, name_cons, rho, young, nu, Gc);
                consMod.push_back(constitutiveModel);
            } else if (line.compare(0, 13, "*ViscoElastic") == 0 && flagNode) {

                // Line 1
                templine.clear();
                getline(ifs, templine);
                double K = atof(templine.c_str());
                // Line 2
                templine.clear();
                getline(ifs, templine);
                int N = atoi(templine.c_str());
                if (N <= 0) {
                    ERROR("The minimum number of branches in this viscoelastic model is 1. You gave %i", N);
                    exit(-1);
                }
                // Line 3. It has several commas
                commas.clear();
                commas.resize(N);
                vector<double> mus(N + 1, 0);
                templine.clear();
                getline(ifs, templine);
                tmp = int(templine.size());
                int contc = 0;
                for (int j = 0; j < tmp; j++) {
                    if (templine.compare(j, 1, ",") == 0) {
                        commas[contc] = j;
                        contc++;
                    }
                }
                templine.copy(tmpme2, commas[0], 0);
                mus[0] = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                for (int j = 1; j < commas.size(); j++) {
                    templine.copy(tmpme2, commas[j] - commas[j - 1] - 1, commas[j - 1] + 1);
                    mus[j] = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                }
                templine.copy(tmpme2, tmp, commas[N - 1] + 1);
                mus[N] = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                // Line 4. It may have several commas
                commas.clear();
                commas.resize(N - 1);
                vector<double> etas(N, 0.);
                templine.clear();
                getline(ifs, templine);
                tmp = int(templine.size());
                if (N == 1) {
                    etas[0] = atof(templine.c_str());
                } else {
                    contc = 0;
                    for (int j = 0; j < tmp; j++) {
                        if (templine.compare(j, 1, ",") == 0) {
                            commas[contc] = j;
                            contc++;
                        }
                    }
                    templine.copy(tmpme2, commas[0], 0);
                    etas[0] = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);

                    for (int j = 1; j < commas.size(); j++) {
                        templine.copy(tmpme2, commas[j] - commas[j - 1] - 1, commas[j - 1] + 1);
                        etas[j] = atof(tmpme2);
                        memset(tmpme2, 0, maxsize);
                    }
                    templine.copy(tmpme2, tmp, commas[N - 1 - 1] + 1);
                    etas[N - 1] = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                }
                name_cons = "ViscoElastic";
                constitutiveModel = new classViscoElastic(ndim, name, name_cons, rho, K, N, mus, etas);
                consMod.push_back(constitutiveModel);
            } else if (line.compare(0, 12, "*ViscoGrowth") == 0 && flagNode) {
                // Line 1
                templine.clear();
                getline(ifs, templine);
                double K = atof(templine.c_str());
                // Line 2
                templine.clear();
                getline(ifs, templine);
                int N = atoi(templine.c_str());
                if (N <= 0) {
                    ERROR("The minimum number of branches in this viscoelastic model is 1. You gave %i", N);
                    exit(-1);
                }
                // Line 3. It has several commas
                commas.clear();
                commas.resize(N);
                vector<double> mus(N + 1, 0);
                templine.clear();
                getline(ifs, templine);
                tmp = int(templine.size());
                int contc = 0;
                for (int j = 0; j < tmp; j++) {
                    if (templine.compare(j, 1, ",") == 0) {
                        commas[contc] = j;
                        contc++;
                    }
                }
                templine.copy(tmpme2, commas[0], 0);
                mus[0] = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                for (int j = 1; j < commas.size(); j++) {
                    templine.copy(tmpme2, commas[j] - commas[j - 1] - 1, commas[j - 1] + 1);
                    mus[j] = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                }
                templine.copy(tmpme2, tmp, commas[N - 1] + 1);
                mus[N] = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                // Line 4. It may have several commas
                commas.clear();
                commas.resize(N - 1);
                vector<double> etas(N, 0.);
                templine.clear();
                getline(ifs, templine);
                tmp = int(templine.size());
                if (N == 1) {
                    etas[0] = atof(templine.c_str());
                } else {
                    contc = 0;
                    for (int j = 0; j < tmp; j++) {
                        if (templine.compare(j, 1, ",") == 0) {
                            commas[contc] = j;
                            contc++;
                        }
                    }
                    templine.copy(tmpme2, commas[0], 0);
                    etas[0] = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);

                    for (int j = 1; j < commas.size(); j++) {
                        templine.copy(tmpme2, commas[j] - commas[j - 1] - 1, commas[j - 1] + 1);
                        etas[j] = atof(tmpme2);
                        memset(tmpme2, 0, maxsize);
                    }
                    templine.copy(tmpme2, tmp, commas[N - 1 - 1] + 1);
                    etas[N - 1] = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                }

                templine.clear();
                getline(ifs, templine);
                double Gc = atof(templine.c_str());
                name_cons = "ViscoGrowth";
                constitutiveModel = new classViscoElasticGrowth(ndim, name, name_cons, rho, K, N, mus, etas, Gc);
                consMod.push_back(constitutiveModel);
            } else if (line.compare(0, 21, "*ViscoMechAxialGrowth") == 0 && flagNode) {
                // Line 1
                templine.clear();
                getline(ifs, templine);
                double K = atof(templine.c_str());
                // Line 2
                templine.clear();
                getline(ifs, templine);
                int N = atoi(templine.c_str());
                if (N <= 0) {
                    ERROR("The minimum number of branches in this viscoelastic model is 1. You gave %i", N);
                    exit(-1);
                }
                // Line 3. It has several commas
                commas.clear();
                commas.resize(N);
                vector<double> mus(N + 1, 0);
                templine.clear();
                getline(ifs, templine);
                tmp = int(templine.size());
                int contc = 0;
                for (int j = 0; j < tmp; j++) {
                    if (templine.compare(j, 1, ",") == 0) {
                        commas[contc] = j;
                        contc++;
                    }
                }
                templine.copy(tmpme2, commas[0], 0);
                mus[0] = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                for (int j = 1; j < commas.size(); j++) {
                    templine.copy(tmpme2, commas[j] - commas[j - 1] - 1, commas[j - 1] + 1);
                    mus[j] = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                }
                templine.copy(tmpme2, tmp, commas[N - 1] + 1);
                mus[N] = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                // Line 4. It may have several commas
                commas.clear();
                commas.resize(N);
                vector<double> etas(N, 0.);
                templine.clear();
                getline(ifs, templine);
                tmp = int(templine.size());
                if (N == 1) {
                    etas[0] = atof(templine.c_str());
                } else {
                    contc = 0;
                    for (int j = 0; j < tmp; j++) {
                        if (templine.compare(j, 1, ",") == 0) {
                            commas[contc] = j;
                            contc++;
                        }
                    }
                    templine.copy(tmpme2, commas[0], 0);
                    etas[0] = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);

                    for (int j = 1; j < commas.size(); j++) {
                        templine.copy(tmpme2, commas[j] - commas[j - 1] - 1, commas[j - 1] + 1);
                        etas[j] = atof(tmpme2);
                        memset(tmpme2, 0, maxsize);
                    }
                    templine.copy(tmpme2, tmp, commas[N - 1 - 1] + 1);
                    etas[N - 1] = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                }

                templine.clear();
                getline(ifs, templine);
                commas.clear();
                commas.resize(3 + ndim);
                vector<double> values(ndim + 4, 0);
                tmp = int(templine.size());
                contc = 0;
                for (int j = 0; j < tmp; j++) {
                    if (templine.compare(j, 1, ",") == 0) {
                        commas[contc] = j;
                        contc++;
                    }
                }
                templine.copy(tmpme2, commas[0], 0);
                values[0] = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                for (int j = 1; j < commas.size(); j++) {
                    templine.copy(tmpme2, commas[j] - commas[j - 1] - 1, commas[j - 1] + 1);
                    values[j] = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                }
                templine.copy(tmpme2, tmp, commas[ndim + 3 - 1] + 1);
                values[ndim + 4 - 1] = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                // Line 4. It may have severa

                double Gc = values[0];// Initial growth multiplier
                double Kp = values[1], Kd = values[2], rhoMicro = values[3];
                vector<double> n0(ndim, 0);
                for (int i = 0; i < ndim; i++) {
                    n0[i] = values[4 + i];
                    //cout<<n0[i]<<endl;
                }
                name_cons = "ViscoMechAxialGrowth";
                constitutiveModel = new classViscoMechAxialGrowth(ndim, name, name_cons, rho, K, N, mus, etas, Gc, n0,
                                                                  Kp, Kd, rhoMicro);
                consMod.push_back(constitutiveModel);
            } else if (line.compare(0, 29, "*ViscoMechContractAxialGrowth") == 0 && flagNode) {
                // Line 1
                templine.clear();
                getline(ifs, templine);
                double K = atof(templine.c_str());
                // Line 2
                templine.clear();
                getline(ifs, templine);
                int N = atoi(templine.c_str());
                if (N <= 0) {
                    ERROR("The minimum number of branches in this viscoelastic model is 1. You gave %i", N);
                    exit(-1);
                }
                // Line 3. It has several commas
                commas.clear();
                commas.resize(N);
                vector<double> mus(N + 1, 0);
                templine.clear();
                getline(ifs, templine);
                tmp = int(templine.size());
                int contc = 0;
                for (int j = 0; j < tmp; j++) {
                    if (templine.compare(j, 1, ",") == 0) {
                        commas[contc] = j;
                        contc++;
                    }
                }
                templine.copy(tmpme2, commas[0], 0);
                mus[0] = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                for (int j = 1; j < commas.size() + 1; j++) {
                    if (j == commas.size()) {
                        templine.copy(tmpme2, templine.size() - commas[j - 1], commas[j - 1] + 1);

                    } else {
                        templine.copy(tmpme2, commas[j] - commas[j - 1] - 1, commas[j - 1] + 1);
                    }
                    mus[j] = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                }
                //	templine.copy(tmpme2,tmp,commas[N-1]+1);
                //	mus[N]=atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                // Line 4. It may have several commas
                commas.clear();
                commas.resize(N);
                vector<double> etas(N, 0.);
                templine.clear();
                getline(ifs, templine);
                tmp = int(templine.size());
                if (N == 1) {
                    etas[0] = atof(templine.c_str());
                } else {
                    contc = 0;
                    for (int j = 0; j < tmp; j++) {
                        if (templine.compare(j, 1, ",") == 0) {
                            commas[contc] = j;
                            contc++;
                        }
                    }
                    templine.copy(tmpme2, commas[0], 0);
                    etas[0] = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);

                    for (int j = 1; j < commas.size(); j++) {
                        templine.copy(tmpme2, commas[j] - commas[j - 1] - 1, commas[j - 1] + 1);
                        etas[j] = atof(tmpme2);
                        memset(tmpme2, 0, maxsize);
                    }
                    templine.copy(tmpme2, tmp, commas[N - 1 - 1] + 1);
                    etas[N - 1] = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                }

                templine.clear();
                getline(ifs, templine);
                commas.clear();
                commas.resize(5 + ndim);
                vector<double> values(ndim + 6, 0);
                tmp = int(templine.size());
                contc = 0;
                for (int j = 0; j < tmp; j++) {
                    if (templine.compare(j, 1, ",") == 0) {
                        commas[contc] = j;
                        contc++;
                    }
                }
                templine.copy(tmpme2, commas[0], 0);
                values[0] = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                for (int j = 1; j < commas.size(); j++) {
                    templine.copy(tmpme2, commas[j] - commas[j - 1] - 1, commas[j - 1] + 1);
                    values[j] = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                }
                templine.copy(tmpme2, tmp, commas[ndim + 5 - 1] + 1);
                values[ndim + 6 - 1] = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                double Gc = values[0];// Initial growth multiplier
                double Kp = values[1], Kd = values[2], rhoMC = values[3], c0 = values[4], S0 = values[5];
                vector<double> n0(ndim, 0);
                for (int i = 0; i < ndim; i++) {
                    n0[i] = values[6 + i];
                }
                name_cons = "ViscoMechContractAxialGrowth";
                constitutiveModel = new classViscoMechContractAxialGrowth(ndim, name, name_cons, rho, K, N, mus, etas,
                                                                          Gc, n0, Kp, Kd, rhoMC, c0, S0);
                consMod.push_back(constitutiveModel);
            } else if (line.compare(0, 13, "*J2Plasticity") == 0 && flagNode) {
                getline(ifs, templine);
                tmp = int(templine.size());
                commas.clear();
                for (int j = 0; j < tmp; j++) {
                    if (templine.compare(j, 1, ",") == 0) {
                        commas.push_back(j);
                    }
                }
                commas.push_back(tmp + 1);
                memset(tmpme2, 0, maxsize);
                templine.copy(tmpme2, commas[0], 0);
                double young = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                templine.copy(tmpme2, commas[1] - commas[0] - 1, commas[0] + 1);
                double nu = atof(tmpme2);
                memset(tmpme2, 0, maxsize);
                templine.copy(tmpme2, commas[2] - commas[1] - 1, commas[1] + 1);
                int order = atoi(tmpme2);
                memset(tmpme2, 0, maxsize);
                commas.clear();
                
                //
                getline(ifs, templine);
                tmp = int(templine.size());
                for (int j = 0; j < tmp; j++) {
                    if (templine.compare(j, 1, ",") == 0) {
                        commas.push_back(j);
                    }
                }
                commas.push_back(tmp + 1);
                memset(tmpme2, 0, maxsize);
                templine.copy(tmpme2, commas[0], 0);
                string hardeningLawName(tmpme2);
                vector<double> data;
                for (int jj=1; jj< commas.size(); jj++)
                {
                    memset(tmpme2, 0, maxsize);
                    templine.copy(tmpme2, commas[jj] - commas[jj-1] - 1, commas[jj-1] + 1);
                    data.push_back(atof(tmpme2));
                }
                IsotropicHardening* hardenLaw = IsotropicHardening::allocate(hardeningLawName.c_str(),data);
                name_cons = "J2Plasticity";
                constitutiveModel = new J2plasticityFiniteStrain(ndim, name, name_cons, rho, young, nu, *hardenLaw, order);
                delete hardenLaw;
                consMod.push_back(constitutiveModel);
            } else if (line.compare(0, 9, "*ExtraDof") == 0 && flagNode) {
                templine.clear();
                getline(ifs, templine);
                if (templine.compare(0, 8, "*FHNelec") == 0) {
                    bool mechaCoupling = false;
                    if (templine.compare(0, 16, "*FHNelecCoupling") == 0) {
                        name_cons = "FHNelecCoupling";
                        mechaCoupling = true;
                    } else {
                        name_cons = "FHNelec";;
                    }

                    getline(ifs, templine);
                    tmp = int(templine.size());
                    commas.clear();
                    int contc = 0;
                    for (int j = 0; j < tmp; j++) {
                        if (templine.compare(j, 1, ",") == 0) {
                            commas.push_back(j);
                            contc++;
                        }
                    }
                    commas.push_back(tmp + 1);
                    templine.copy(tmpme2, commas[0], 0);
                    vector<double> n0(3, 0);
                    n0[0] = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                    templine.copy(tempNode, commas[1] - commas[0] - 1, commas[0] + 1);
                    n0[1] = atof(tempNode);
                    memset(tempNode, 0, maxsize);
                    templine.copy(tempNode, commas[2] - commas[1] - 1, commas[1] + 1);
                    n0[2] = atof(tempNode);
                    memset(tempNode, 0, maxsize);
                    templine.copy(tempNode, commas[3] - commas[2] - 1, commas[2] + 1);
                    double D_iso = 0;
                    D_iso = atof(tempNode);
                    memset(tempNode, 0, maxsize);
                    templine.copy(tempNode, commas[4] - commas[3] - 1, commas[3] + 1);
                    double D_ani = 0;
                    D_ani = atof(tempNode);
                    memset(tempNode, 0, maxsize);
                    templine.copy(tempNode, commas[5] - commas[4] - 1, commas[4] + 1);
                    double aelec = 0;
                    aelec = atof(tempNode);
                    memset(tempNode, 0, maxsize);
                    templine.copy(tempNode, commas[6] - commas[5] - 1, commas[5] + 1);
                    double b = 0;
                    b = atof(tempNode);
                    memset(tempNode, 0, maxsize);
                    templine.copy(tempNode, commas[7] - commas[6] - 1, commas[6] + 1);
                    double epsi = 0;
                    epsi = atof(tempNode);
                    memset(tempNode, 0, maxsize);
                    templine.copy(tempNode, commas[8] - commas[7] - 1, commas[7] + 1);
                    double vequi = 0;
                    vequi = atof(tempNode);
                    memset(tempNode, 0, maxsize);
                    templine.copy(tempNode, commas[9] - commas[8] - 1, commas[8] + 1);
                    double ionicCurrent = 0;
                    ionicCurrent = atof(tempNode);
                    memset(tempNode, 0, maxsize);
                    commas.clear();
                    constitutiveModel = new classFHN(ndim, name, name_cons, rho, aelec, b, epsi, vequi, n0, D_iso,
                                                     D_ani, ionicCurrent, mechaCoupling);
                    consMod.push_back(constitutiveModel);
                } else if (templine.compare(0, 12, "*Temperature") == 0) {
                    bool mechaCoupling=false;
                    if (templine.compare(0, 20, "*TemperatureCoupling") == 0) {
                        name_cons = "CouplingTemperature";
                        mechaCoupling = true;
                    } else {
                        name_cons = "Temperature";
                    }
                    getline(ifs, templine);
                    tmp = int(templine.size());
                    commas.clear();
                    int contc = 0;
                    for (int j = 0; j < tmp; j++) {
                        if (templine.compare(j, 1, ",") == 0) {
                            commas.push_back(j);
                            contc++;
                        }
                    }
                    commas.push_back(tmp + 1);
                    templine.copy(tmpme2, commas[0], 0);
                    double Cc_0;
                    Cc_0 = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                    templine.copy(tempNode, commas[1] - commas[0] - 1, commas[0] + 1);
                    double Cc_1;
                    Cc_1 = atof(tempNode);
                    memset(tempNode, 0, maxsize);
                    templine.copy(tempNode, commas[2] - commas[1] - 1, commas[1] + 1);
                    double k_0;
                    k_0 = atof(tempNode);
                    memset(tempNode, 0, maxsize);
                    templine.copy(tempNode, commas[3] - commas[2] - 1, commas[2] + 1);
                    double k_1;
                    k_1 = atof(tempNode);
                    memset(tempNode, 0, maxsize);
                    templine.copy(tempNode, commas[4] - commas[3] - 1, commas[3] + 1);
                    double theta_ini;
                    theta_ini = atof(tempNode);
                    memset(tempNode, 0, maxsize);
                    templine.copy(tempNode, commas[5] - commas[4] - 1, commas[4] + 1);
                    double theta_ref1;
                    theta_ref1 = atof(tempNode);
                    memset(tempNode, 0, maxsize);
                    templine.copy(tempNode, commas[6] - commas[5] - 1, commas[5] + 1);
                    double theta_ref2;
                    theta_ref2 = atof(tempNode);
                    memset(tempNode, 0, maxsize);
                    commas.clear();
                    constitutiveModel = new classtemperature(ndim, name, name_cons, rho, Cc_0, Cc_1, k_0, k_1,
                                                             theta_ini, theta_ref1, theta_ref2, mechaCoupling);
                    consMod.push_back(constitutiveModel);
                } else if (templine.compare(0, 22, "*3DprintingTemperature") == 0) {
                    bool mechaCoupling=false;
                    if (templine.compare(0, 30, "*3DprintingTemperatureCoupling") == 0) {
                        name_cons = "3DprintingTemperatureCoupling";
                        mechaCoupling = true;
                    } else {
                        name_cons = "3DprintingTemperature";
                    }
                    getline(ifs, templine);
                    tmp = int(templine.size());
                    commas.clear();
                    int contc = 0;
                    for (int j = 0; j < tmp; j++) {
                        if (templine.compare(j, 1, ",") == 0) {
                            commas.push_back(j);
                            contc++;
                        }
                    }
                    commas.push_back(tmp + 1);
                    templine.copy(tmpme2, commas[0], 0);
                    double CRefSolid = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                    templine.copy(tmpme2, commas[1] - commas[0] - 1, commas[0] + 1);
                    double C1_SL = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                    templine.copy(tmpme2, commas[2] - commas[1] - 1, commas[1] + 1);
                    double C2_SL = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                    templine.copy(tmpme2, commas[3] - commas[2] - 1, commas[2] + 1);
                    double latentHeat1_PL = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                    templine.copy(tmpme2, commas[4] - commas[3] - 1, commas[3] + 1);
                    double latentHeat2_PL = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                    templine.copy(tmpme2, commas[5] - commas[4] - 1, commas[4] + 1);
                    double latentHeat3_PL = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                    templine.copy(tmpme2, commas[6] - commas[5] - 1, commas[5] + 1);
                    double latentHeat1_SL = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                    templine.copy(tmpme2, commas[7] - commas[6] - 1, commas[6] + 1);
                    double latentHeat2_SL = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                    templine.copy(tmpme2, commas[8] - commas[7] - 1, commas[7] + 1);
                    double latentHeat3_SL = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                    templine.copy(tmpme2, commas[9] - commas[8] - 1, commas[8] + 1);
                    double kRefPowder = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
		            templine.copy(tmpme2, commas[10] - commas[9] - 1, commas[9] + 1);
                    double k1_PL = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                    templine.copy(tmpme2, commas[11] - commas[10] - 1, commas[10] + 1);
                    double kRefSolid = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                    templine.copy(tmpme2, commas[12] - commas[11] - 1, commas[11] + 1);
                    double k1_SL = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                    templine.copy(tmpme2, commas[13] - commas[12] - 1, commas[12] + 1);
                    double k2_SL = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
		            templine.copy(tmpme2, commas[14] - commas[13] - 1, commas[13] + 1);
                    double k3_SL = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                    templine.copy(tmpme2, commas[15] - commas[14] - 1, commas[14] + 1);
                    double logGrowthRate = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                    templine.copy(tmpme2, commas[16] - commas[15] - 1, commas[15] + 1);
                    double  theta_sigMidpoint = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
		            templine.copy(tmpme2, commas[17] - commas[16] - 1, commas[16] + 1);
                    double theta_ref = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                    templine.copy(tmpme2, commas[18] - commas[17] - 1, commas[17] + 1);
                    double theta_powder = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                    templine.copy(tmpme2, commas[19] - commas[18] - 1, commas[18] + 1);
                    double theta_SL = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                    templine.copy(tmpme2, commas[20] - commas[19] - 1, commas[19] + 1);
                    double theta_liq = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                    commas.clear();
		            constitutiveModel = new classtemperature3Dprinting(ndim, name, name_cons, rho, CRefSolid,
                                                                       C1_SL, C2_SL, latentHeat1_PL,
                                                                       latentHeat2_PL,
                                                                       latentHeat3_PL, latentHeat1_SL, latentHeat2_SL,
                                                                       latentHeat3_SL, kRefPowder, k1_PL,
                                                                       kRefSolid, k1_SL, k2_SL, k3_SL,
								       logGrowthRate, theta_sigMidpoint,
                                                                       theta_ref, theta_powder, theta_SL,
                                                                       theta_liq, mechaCoupling);
                    consMod.push_back(constitutiveModel);
                }
                flagExtra = true;
                // Getting a0 (anisotropic direction) of the elements
                templine.clear();
                streampos oldpos;
                oldpos = ifs.tellg();
                while (getline(ifs, templine)) {
                    if (templine.compare(0, 1, "*") == 0) {
                        ifs.seekg(oldpos);
                        break;
                    } else {
                        int sizeLine = int(templine.size());
                        commas.clear();
                        int contc = 0;
                        for (int j = 0; j < sizeLine; j++) {
                            if (templine.compare(j, 1, ",") == 0) {
                                commas.push_back(j);
                                contc++;
                            }
                        }
                        commas.push_back(sizeLine + 1);

                        int elem = 0;
                        vector<double> a0(3, 0);

                        memset(tmpme2, 0, maxsize);
                        templine.copy(tmpme2, commas[0], 0);
                        elem = atoi(tmpme2) - 1;
#ifdef PARALLEL
                        if (epart[elem] == rank) {
                            it = find(elemLocal.begin(), elemLocal.end(), elem);
                            elem = it - elemLocal.begin();

                            memset(tmpme2, 0, maxsize);
                            templine.copy(tmpme2, commas[1] - commas[0] - 1, commas[0] + 1);
                            a0[0] = atof(tmpme2);
                            memset(tmpme2, 0, maxsize);
                            templine.copy(tmpme2, commas[2] - commas[1] - 1, commas[1] + 1);
                            a0[1] = atof(tmpme2);
                            memset(tmpme2, 0, maxsize);
                            templine.copy(tmpme2, commas[3] - commas[2] - 1, commas[2] + 1);
                            a0[2] = atof(tmpme2);
                            memset(tmpme2, 0, maxsize);
                            elements[elem]->set_n0(a0);
                        }
#else
                        memset(tmpme2, 0, maxsize);
                        templine.copy(tmpme2, commas[1] - commas[0] - 1, commas[0] + 1);
                        a0[0] = atof(tmpme2);
                        memset(tmpme2, 0, maxsize);
                        templine.copy(tmpme2, commas[2] - commas[1] - 1, commas[1] + 1);
                        a0[1] = atof(tmpme2);
                        memset(tmpme2, 0, maxsize);
                        templine.copy(tmpme2, commas[3] - commas[2] - 1, commas[2] + 1);
                        a0[2] = atof(tmpme2);
                        memset(tmpme2, 0, maxsize);
                        elements[elem]->set_n0(a0);
#endif
                    }
                    oldpos = ifs.tellg();
                }
            } else if (line.compare(0, 11, "*Stochastic") == 0 && flagNode) {
                flagStochastic = true;
                templine.clear();
                getline(ifs, templine);
                if (templine.compare(0, 18, "*StochasticMapping") == 0) {                  
                    stochasticMapping = true;
                    templine.clear();
                    getline(ifs, templine);
                }
                if (templine.compare(0, 13, "*HyperElastic") == 0 && flagNode) {
                    if ((!templine.compare(15, 34, "St-Venant-Kirchhoff") == 0) && (!templine.compare(15, 26, "Neo-Hookean") == 0) &&
                    (!templine.compare(15, 42, "3D-Compressible-Neo-Hookean") == 0)) {
                        ERROR("The %s model is not implemented in MuPhiSim", line.c_str());
                        exit(-1);
                    }
                    string stochModel = templine;
                    printf("%s\n",stochModel.c_str());
                    getline(ifs, templine);
                    tmp = int(templine.size());
                    int contc = 0;
                    for (int j = 0; j < tmp; j++) {
                        if (templine.compare(j, 1, ",") == 0) {
                            commas.push_back(j);
                            contc++;
                        }
                    }
                    commas.push_back(tmp + 1);
                    templine.copy(tmpme2, commas[0], 0);
                    double young = atof(tmpme2);
                    memset(tmpme2, 0, maxsize);
                    templine.copy(tempNode, commas[1] - commas[0] - 1, commas[0] + 1);
                    double nu = atof(tempNode);
                    memset(tempNode, 0, maxsize);
                    int stressState = 0;//default plane stress
                    if (commas.size() == 3)// non default
                    {
                        templine.copy(tempNode, commas[2] - commas[1] - 1, commas[1] + 1);
                        string stressStateString = "";

                        for (int i = 0; i < sizeof(tempNode) / sizeof(char); i++) {
                            if (tempNode[i] != ' ')//remove all whitespaces
                            {
                                stressStateString = stressStateString + tempNode[i];
                            }
                        }
                        if (stressStateString.compare(0, 11, "PlaneStress") == 0) {
                            stressState = 0;
                        } else if (stressStateString.compare(0, 11, "PlaneStrain") == 0) {
                            stressState = 1;
                        } else {
                            stressState = 0;
                        }
                        memset(tempNode, 0, maxsize);
                    }

                    commas.clear();
                    getline(ifs, templine);
                    int order;
                    int resolution = 0;
                    if (templine.compare(0, 11, "*RESOLUTION") == 0) {
                        getline(ifs, templine);
                        templine.copy(tempNode, templine.size(), 0);
                        resolution = atof(tempNode);
                        getline(ifs, templine);
                    }
                    if (templine.compare(0, 6, "*ORDER") == 0) {
                        getline(ifs, templine);
                        templine.copy(tempNode, templine.size(), 0);
                        order = atof(tempNode);
                    }
                    getline(ifs, templine);
                    vector<double> stochastic_function;
                    vector<double> stochastic_parameters;
                    vector<double> stochastic_number;
                    int cont_stochastic_function;
                    while (!(templine.compare(0, 11, "*POLYNOMIAL") == 0)) {
                        cont_stochastic_function = 0;
                        tmp = int(templine.size());
                        contc = 0;
                        for (int j = 0; j < tmp; j++) {
                            if (templine.compare(j, 1, ",") == 0) {
                                commas.push_back(j);
                                contc++;
                            }
                        }
                        commas.push_back(tmp + 1);
                        templine.copy(tmpme2, commas[0], 0);
                        stochastic_parameters.push_back(atof(tmpme2));
                        memset(tmpme2, 0, maxsize);
                        templine.copy(tmpme2, commas[1] - commas[0] - 1, commas[0] + 1);
                        stochastic_number.push_back(atof(tmpme2));
                        memset(tmpme2, 0, maxsize);
                        for (int i = 0; i < stochastic_number[cont_stochastic_function]; i++) {
                            templine.copy(tempNode, commas[2 + i] - commas[1 + i] - 1, commas[1 + i] + 1);
                            stochastic_function.push_back(atof(tempNode));
                            memset(tempNode, 0, maxsize);
                        }
                        cont_stochastic_function++;
                        getline(ifs, templine);
                        commas.clear();
                    }
                    getline(ifs, templine);
                    string approximation = "";
                    for (int i = 0; i < templine.size(); i++) {
                        if (templine[i] != ' ')//remove all whitespaces
                        {
                            approximation = approximation + templine[i];
                        }
                    }
                    getline(ifs, templine);
                    getline(ifs, templine);
                    string distribution = "";
                    for (int i = 0; i < templine.size(); i++) {
                        if (templine[i] != ' ')//remove all whitespaces
                        {
                            distribution = distribution + templine[i];
                        }
                    }
                    if((!distribution.compare(0, 7, "UNIFORM") == 0) &&
                    (!distribution.compare(0, 8, "GAUSSIAN") == 0)) {
                        WARNING("Only UNIFORM and GAUSSIAN distributions are taken into account at the moment, by default we suppose in this computation the distribution is uniform");
                        distribution = "UNIFORM";
                    }
                    printf("%s\n",stochModel.c_str());
                    if (stochModel.compare(15, 34, "St-Venant-Kirchhoff") == 0) {
                    name_cons = "Stochastic-St-Venant-Kirchhoff";
                    INFO("The stochastic constitutive model is St-Venant-Kirchoff");
                    constitutiveModel = new classStochasticHyperElasticStVenantKirchhoff(ndim, name, name_cons, rho,
                                                                                         young, nu, stressState,
                                                                                         stochastic_number,
                                                                                         stochastic_parameters,
                                                                                         stochastic_function,
                                                                                         approximation, order,
                                                                                         resolution, stochasticMapping, distribution);
                    } else if (stochModel.compare(15, 26, "Neo-Hookean") == 0) {
                    name_cons = "NeoHookean";
                    INFO("The stochastic constitutive model is Neo-Hookean");
                    constitutiveModel = new classStochasticHyperElasticNeoHookean(ndim, name, name_cons, rho,
                                                                                         young, nu,
                                                                                         stochastic_number,
                                                                                         stochastic_parameters,
                                                                                         stochastic_function,
                                                                                         approximation, order,
                                                                                         resolution, stochasticMapping, distribution);    
                    } else if (stochModel.compare(15, 42, "3D-Compressible-Neo-Hookean") == 0) {
                    name_cons = "3D-Compressible-Neo-Hookean";
                    INFO("The stochastic constitutive model is 3D-Compressible-Neo-Hookean");
                    constitutiveModel = new classStochasticHyperElastic3DCompressibleNeoHookean(ndim, name, name_cons, rho,
                                                                                         young, nu,
                                                                                         stochastic_number,
                                                                                         stochastic_parameters,
                                                                                         stochastic_function,
                                                                                         approximation, order,
                                                                                         resolution, stochasticMapping, distribution);    
                    }
                    consMod.push_back(constitutiveModel);
                }

            } else if (line.compare(0, 17, "*List of Elements") == 0 && flagNode) {
                getline(ifs, templine);
                if (constitutiveModel == NULL) {
                    ERROR("The constitutive model has not been properly defined");
                    WARNING("MuPhiSim is a case sensitive program");
                    exit(-1);
                }

                if (templine.compare(0, 3, "All") == 0) {
                    int nEle = elements.size();
                    for (int i = 0; i < nEle; i++) {
                        elements[i]->setFlagExtra(flagExtra);
                        elements[i]->setFlagStochastic(flagStochastic);
                        elements[i]->setConsModel(consMod);


                    }
                } else {
                    flagList = true;
                    size = int(templine.size());
                    templine.copy(tempNode, size, 0);
                    tmp = atoi(tempNode) - 1;
                    listElements.push_back(
                            tmp);//We create a list with the elements in which we have assigned the material
                    if (epart.size() == 0) {
                        elements[tmp]->setFlagExtra(flagExtra);
                        elements[tmp]->setFlagStochastic(flagStochastic);
                        elements[tmp]->setConsModel(consMod);

                    } else {
                        if (epart[tmp] == rank) {
                            it = find(elemLocal.begin(), elemLocal.end(), tmp);
                            tmp = it - elemLocal.begin();
                            elements[tmp]->setFlagExtra(flagExtra);
                            elements[tmp]->setFlagStochastic(flagStochastic);
                            if (flagExtra) {
                                elements[tmp]->setConsModel(consMod);

                            } else {
                                elements[tmp]->setConsModel(constitutiveModel);
                            }
                        }
                    }
                    memset(tempNode, 0, maxsize);
                }
            } else if (line.compare(0, 12, "*INITINTVARS") == 0) {
                //These lines allow to define the value of the internal variables used in the chosen constitutive model
                templine.clear();
                streampos oldpos;
                //Next loop allows to read all the lines that define the internal variables until reach a new keyword (*)
                while (getline(ifs, templine)) {

                    if (templine.compare(0, 1, "*") == 0) { //the loop ends if it has reached * a new keyword
                        ifs.seekg(
                                oldpos); //the code moves to the previous line before reaching "*" (that's necessary to read the keyword again in the global loop)
                        break;
                    } else { //if not *, it reads the information to define the value of a specific internal variable
                        tmp = int(templine.size());
                        commas.clear();
                        int contc = 0;
                        for (int j = 0; j < tmp; j++) {
                            if (templine.compare(j, 1, ",") == 0) {
                                commas.push_back(j);
                                contc++;
                            }
                        }
                        commas.push_back(tmp + 1);
                        memset(tmpme2, 0,
                               maxsize); //Without it, there is a problem if element 1 is not activated at the beginning
                        templine.copy(tmpme2, commas[0], 0);
                        int elemAct = (atof(tmpme2) -
                                       1); //Element in which the internal variable is defined. In c ++ the first element should be zero
                        memset(tmpme2, 0, maxsize);
                        if (flagList) {
                            if (std::find(listElements.begin(), listElements.end(), elemAct) ==
                                listElements.end()) {//Check if the element is in this *MATERIAL Region
                                ERROR("Element on INITINTVARS list not valid");
                                exit(-1);
                            }
                        };
                        templine.copy(tmpme2, commas[1] - commas[0] - 1, commas[0] + 1);
                        int intVarPosition = (atof(tmpme2) -
                                              1); //Internal variable position. In c++ the first position is the zero
                        memset(tmpme2, 0, maxsize);
                        templine.copy(tmpme2, commas[2] - commas[1] - 1, commas[1] + 1);
                        double intVarValue = atof(tmpme2); //Internal variable value
                        memset(tmpme2, 0, maxsize);
                        commas.clear();
                        if (epart.size() == 0) { //if the code is running in sequential
                            elements[elemAct]->add_IntVarsPosition(intVarPosition);
                            elements[elemAct]->add_IntVarsValue(intVarValue);

                        } else {
                            if (epart[elemAct] == rank) { //if the code is running in parallel
                                it = find(elemLocal.begin(), elemLocal.end(), elemAct);
                                elemAct = it - elemLocal.begin();
                                elements[elemAct]->add_IntVarsPosition(intVarPosition);
                                elements[elemAct]->add_IntVarsValue(intVarValue);
                            }
                        }
                    }
                    oldpos = ifs.tellg(); //previous position is stored before get a new line
                }

            } else if (line.compare(0, 13, "*Bounding box") == 0 && flagNode) {
                getline(ifs, templine);
                if (constitutiveModel == NULL) {
                    ERROR("The constitutive model has not been properly defined");
                    WARNING("MuPhiSim is a case sensitive program");
                    exit(-1);
                }
                ERROR("The bounding box option for the constitutive models distribution is not ready yet");
                exit(-1);

                if (templine.compare(0, 3, "All") == 0) {
                    int nEle = elements.size();
                    for (int i = 0; i < nEle; i++) {
                        elements[i]->setConsModel(constitutiveModel);
                    }
                } 
            } else if ((line.compare(0, 4, "*END") == 0 || line.compare(0, 4, "*End") == 0) && flagNode) {
                break;
            } else {
                if (flagList) {
                    if (constitutiveModel == NULL) {
                        ERROR("The constitutive model has not been properly defined");
                        WARNING("MuPhiSim is a case sensitive program");
                        exit(-1);
                    }
                    size = int(line.size());
                    line.copy(tempNode, size, 0);
                    tmp = atoi(tempNode) - 1;
                    listElements.push_back(
                            tmp);//We create a list with the elements in which we have assigned the material
                    if (epart.size() == 0) {
                        elements[tmp]->setFlagExtra(flagExtra);
                        elements[tmp]->setFlagStochastic(flagStochastic);

                        if (flagExtra) {
                            elements[tmp]->setConsModel(consMod);

                        } else {
                            elements[tmp]->setConsModel(constitutiveModel);
                        }
                    } else {
                        if (epart[tmp] == rank) {
                            it = find(elemLocal.begin(), elemLocal.end(), tmp);
                            tmp = it - elemLocal.begin();
                            elements[tmp]->setFlagExtra(flagExtra);
                            elements[tmp]->setFlagStochastic(flagStochastic);

                            if (flagExtra) {
                                elements[tmp]->setConsModel(consMod);

                            } else {
                                elements[tmp]->setConsModel(constitutiveModel);
                            }
                        }
                    }
                    memset(tempNode, 0, maxsize);
                }
            }
        }
    }
    if (!flagNode) {
        ERROR("There is not material models in the input file");
        exit(-1);
    }
}


/*! \brief It assigns the initial value to extraDof variable (it is an optional option. If it is not defined in the input file, extradofs are initialised to zero)
  @param[in] ifs Pointer to the input file (.inp)
  @param[out] iniValueExtraDof Initial value for extradof
*/
extern void initialConditionsManagement(ifstream &ifs, vector<double> &iniValueExtraDof) {

    int maxsize = 100;
    char tempNode[maxsize];
    memset(tempNode, 0, maxsize);
    string line, templine;
    char tmpme2[maxsize];
    memset(tmpme2, 0, maxsize);
    int tmp = 0;
    vector<int> commas;


    while (getline(ifs, line)) {
        if (line.compare(0, 28, "*INITIAL CONDITIONS EXTRADOF") == 0) {
            getline(ifs, templine);
            tmp = int(templine.size());
            int contc = 0;
            for (int j = 0; j < tmp; j++) {
                if (templine.compare(j, 1, ",") == 0) {
                    commas.push_back(j);
                    contc++;
                }
            }
            commas.push_back(tmp + 1);

            for (int j = 0; j < commas.size(); j++) {
                memset(tmpme2, 0, maxsize);
                if (j == 0) {
                    templine.copy(tmpme2, commas[j], 0); //first vale
                } else {
                    templine.copy(tmpme2, commas[j + 1] - commas[j - 1] - 1, commas[j - 1] + 1); //consecutive values
                }
                iniValueExtraDof.push_back(atof(tmpme2)); //vector that stores the initial value of extradof
            }
            break;
        }
    }
#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD); //Wait all the processors to reach this point
#endif
}

extern void dataExtractionManagement(ifstream &ifs, vector<classElements *> &elements,
                             vector<classNodes *> &nodes, vector<int> &npart,
                             vector<int> &nodLocal, vector<int> &epart, vector<int> &elemLocal)
{
    string line;    
    while (getline(ifs, line)) {
        if (line.compare(0, 17, "*EXTRACTION, NODE") == 0) {
            // data extraction
            NodeDataExtractionManagement(nodes,ifs,npart,nodLocal);
        } else if (line.compare(0, 20, "*EXTRACTION, ELEMENT") == 0) {
            // data extraction
            ElementDataExtractionManagement(elements,ifs,epart,elemLocal);
        }
    }
}

/*! \brief It reads the solver and the scale factor.
  @param[out] solver Array with all solvers to be executed
  @param[in] ifs Pointer to the input file (.inp)
  @param[in] elements Array with all elements (or background cells) in the domain
  @param[in] nodes Array with all nodes in the domain
  @param[in] nDof number of dof
  @param[in] ndim dimension of the domain
  @param[in] npart Processor ownership mapping of the nodes in the whole model npart[node i]=rank (empty for sequential simulations)
  @param[in] nodLocal Mapping of local nodes to global nodes numbering nodlocal[local node index]= global node index (empty for sequential simulations)
*/
extern void solverManagement(vector<solvers *> &solver, ifstream &ifs, vector<classElements *> &elements,
                             vector<classNodes *> &nodes, int nDof, int ndim, vector<int> &npart,
                             vector<int> &nodLocal, vector<int> &epart, vector<int> &elemLocal) {

    streampos oldpos;
    int maxsize = 100;
    char tempNode[maxsize];
    string name;
    memset(tempNode, 0, maxsize);
    bool flagSolver = false;
    bool flagFSI = false;
    string line, line2, secondLine, templine, templine2;
    int size = 0;
    char tmpme2[maxsize];
    memset(tmpme2, 0, maxsize);
    vector<int> commas;
    string solverType;
    double scaleFactor = 1;
    int numSteps = 0;
    double relTol=-1.;
    double absTol=-1.;
    int maximalNumberIterations = 0;
    vector<classDirichlet *> DNod; /* Dirichlet (essential) boundary conditions*/
    double t0 = 0, tf = 0;
    classOutputs *outputs=NULL;
    classPrintForces *forceOutputs = nullptr;
    solvers *sol=NULL;
    vector<int> activeDofMapLocal;
    vector<int> nActiveDof;
    vector<int> ISactiveDof;
    bool tangentByPerturbation = false;
    double perturbationTol = 1e-8;
    TimeStepping steppingPlan(1e-12);
    string timeStepAdaptiveType="";
    vector<double> timeStepAdaptiveData;
    double maximalTimeStep = -1.;
    /// Petsc solver variables
    string petsc_precon = "none";
    string petsc_solver = "lu";
    int flag_solver = 0;
    int flag_pc = 0;
    int allowableNumFailedSteps=0;
    
    string trueDispField ="";
    vector<double> trueDispFieldData;
    
    solvers::StiffnessMatrixType stiffnessMatrixType=solvers::None;
    int fixedMatrixIteration = -1;

    vector<double> val_viscosity;
    val_viscosity.resize(2);
    val_viscosity[0] = 0.;
    val_viscosity[1] = 0.;
    int flag_AV = 0;
    bool flagActivationSolver = false;  //This flag active the option of birth/death elements
    bool flagInstantaneous = false;
    //flagInstantaneous = true Dirichlet BC are applied instantaneously, flagInstantaneous = false if they are applied as a ramp

    while (getline(ifs, line)) {
        if (line.compare(0, 7, "*SOLVER") == 0) {
            flagSolver = true;
            if (line.compare(9, 18, "EXPLICIT") == 0) {
                solverType = "EXPLICIT";
                val_viscosity[0] = 0.;
                val_viscosity[1] = 0.;
                flag_AV = 0;
            } else if (line.compare(9, 18, "IMPLICIT") == 0) {
                solverType = "IMPLICIT";
            } else if (line.compare(9, 25, "IMPLICIT STATIC") == 0) {
                solverType = "IMPLICIT STATIC";
            } else {
                ERROR("The solver has not been properly defined");
                WARNING("MuPhiSim is a case sensitive program");
                exit(-1);
            }
            while (getline(ifs, line2)) {
                if (line2.compare(0, 12, "*END SOLVER") == 0) {
                    break;
                } else {
                    if (line2.compare(0, 11, "*ACTIVATION") == 0) {
                        //This line actives the option of birth/death elements
                        flagActivationSolver = true;
                    } else if (line2.compare(0, 4, "*FSI") == 0) {
                        flagFSI = true;
                        INFO("FSI simulation");
                    } else if (line2.compare(0, 13, "*Scale Factor") == 0) {
                        templine2.clear();
                        getline(ifs, templine2);
                        size = int(templine2.size());
                        templine2.copy(tempNode, size, 0);
                        scaleFactor = atof(tempNode);
                        memset(tempNode, 0, maxsize);
                        templine2.clear();
                    } else if (line2.compare(0, 16, "*Number Of Steps") == 0) {
                        templine2.clear();
                        getline(ifs, templine2);
                        size = int(templine2.size());
                        templine2.copy(tempNode, size, 0);
                        numSteps = atoi(tempNode);
                        memset(tempNode, 0, maxsize);
                        templine2.clear();
                        steppingPlan.setNumberOfSteps(numSteps);
                     } else if (line2.compare(0, 33, "*Allowable Number of Failed Steps") == 0) {
                        templine2.clear();
                        getline(ifs, templine2);
                        size = int(templine2.size());
                        templine2.copy(tempNode, size, 0);
                        allowableNumFailedSteps = atoi(tempNode);
                        memset(tempNode, 0, maxsize);
                        templine2.clear();
                     } else if (line2.compare(0, 14, "*Max Time Step") == 0) {
                        templine2.clear();
                        getline(ifs, templine2);
                        size = int(templine2.size());
                        templine2.copy(tempNode, size, 0);
                        maximalTimeStep = atof(tempNode);
                        memset(tempNode, 0, maxsize);
                        templine2.clear();
                     } else if (line2.compare(0, 19, "*Absolute Tolerance") == 0) {
                        templine2.clear();
                        getline(ifs, templine2);
                        size = int(templine2.size());
                        templine2.copy(tempNode, size, 0);
                        absTol = atof(tempNode);
                        memset(tempNode, 0, maxsize);
                        templine2.clear();
                    } else if (line2.compare(0, 19, "*Relative Tolerance") == 0) {
                        templine2.clear();
                        getline(ifs, templine2);
                        size = int(templine2.size());
                        templine2.copy(tempNode, size, 0);
                        relTol = atof(tempNode);
                        memset(tempNode, 0, maxsize);
                        templine2.clear();
                    } else if (line2.compare(0, 15, "*Max Iterations") == 0) {
                        templine2.clear();
                        getline(ifs, templine2);
                        size = int(templine2.size());
                        templine2.copy(tempNode, size, 0);
                        maximalNumberIterations = atof(tempNode);
                        memset(tempNode, 0, maxsize);
                        templine2.clear();
                    } else if (line2.compare(0, 24, "*TANGENT BY PERTURBATION") == 0) {
                        tangentByPerturbation = true;
                        templine2.clear();
                        getline(ifs, templine2);
                        size = int(templine2.size());
                        templine2.copy(tempNode, size, 0);
                        perturbationTol = atof(tempNode);
                        memset(tempNode, 0, maxsize);
                        templine2.clear();
                    } else if (line2.compare(0, 15, "*ERROR ANALYSIS") == 0) {
                        commas.clear();
                        //
                        templine.clear();
                        getline(ifs, templine);
                        size = int(templine.size());
                        for (int j = 0; j < size; j++) {
                            if (templine.compare(j, 1, ",") == 0) {
                                commas.push_back(j);
                            }
                        }
                        commas.push_back(size + 1);
                        memset(tempNode, 0, maxsize);
                        templine.copy(tempNode, commas[0], 0);
                        trueDispField = tempNode;
                        for (int jj=1; jj< commas.size(); jj++)
                        {
                            memset(tempNode, 0, maxsize);
                            templine.copy(tempNode, commas[jj] - commas[jj-1] - 1, commas[jj-1] + 1);
                            trueDispFieldData.push_back(atof(tempNode));
                        }                    
                    } else if (line2.compare(0, 17, "*STIFFNESS MATRIX") == 0) {
                        int tmp2 = int(line2.size());
                        if (line2.compare(24, tmp2, "CURRENT") == 0) 
                        {
                            stiffnessMatrixType = solvers::Current;
                        }
                        else if (line2.compare(24, tmp2, "INITIAL") == 0)  
                        {
                            stiffnessMatrixType = solvers::Initial;
                        }
                        else if (line2.compare(24, tmp2, "ITERATION") == 0)  
                        {
                            stiffnessMatrixType = solvers::Iteration;
                            templine2.clear();
                            getline(ifs, templine2);
                            size = int(templine2.size());
                            templine2.copy(tempNode, size, 0);
                            fixedMatrixIteration = atof(tempNode);
                            memset(tempNode, 0, maxsize);
                            templine2.clear();
                        }
                        else if (line2.compare(24, tmp2, "INTERVAL") == 0)  
                        {
                            stiffnessMatrixType = solvers::Interval;
                            templine2.clear();
                            getline(ifs, templine2);
                            size = int(templine2.size());
                            templine2.copy(tempNode, size, 0);
                            fixedMatrixIteration = atof(tempNode);
                            memset(tempNode, 0, maxsize);
                            templine2.clear();
                        }
                        else if (line2.compare(24, tmp2, "BEGINNING") == 0)  
                        {
                            stiffnessMatrixType = solvers::Beginning;
                        } 
         
                    } else if (line2.compare(0, 15, "*Bulk viscosity") == 0) {
                        val_viscosity.resize(2);
                        val_viscosity[0] = 0.06;
                        val_viscosity[1] = 1.2;
                        string tempChar;
                        flag_AV = 1;
                        templine.clear();
                        getline(ifs, templine);
                        int tmp = int(templine.size());
                        for (int i = 0; i < tmp; i++) {
                            if (templine.compare(i, 1, "*") == 0) { //No parameter
                                break;
                            } else if (templine.compare(i, 1, " ") != 0) {
                                if ((templine.compare(i, 1, ",") == 0)) {
                                    val_viscosity[0] = atof(tempChar.c_str());
                                    tempChar.clear();
                                } else {
                                    tempChar += templine[i];
                                    if (i == tmp - 1) {
                                        val_viscosity[1] = atof(tempChar.c_str());
                                        tempChar.clear();
                                    }
                                }
                            }
                        }
                    } else if (line2.compare(0, 14, "*Dev viscosity") == 0) {
                        val_viscosity.resize(2);
                        val_viscosity[0] = 1;
                        val_viscosity[1] = 0.5;
                        string tempChar = "";
                        flag_AV = 2;
                        templine.clear();
                        getline(ifs, templine);
                        int tmp = int(templine.size());
                        for (int i = 0; i < tmp; i++) {
                            if (templine.compare(i, 1, "*") == 0) { //No parameter
                                break;
                            } else if (templine.compare(i, 1, " ") != 0) {
                                if ((templine.compare(i, 1, ",") == 0)) {
                                    val_viscosity[0] = atof(tempChar.c_str());
                                    tempChar.clear();
                                } else {
                                    tempChar += templine[i];
                                    if (i == tmp - 1) {
                                        val_viscosity[1] = atof(tempChar.c_str());
                                        tempChar.clear();
                                    }
                                }
                            }
                        }
                    } else if (line2.compare(0, 33, "*BOUNDARY, TYPE=DISPLACEMENT RAMP") == 0) {

#ifdef PARALLEL
                        MPI_Barrier(MPI_COMM_WORLD); //Wait all the processors to reach this point
#endif
                        double duration = tf - t0;
                        vector<constitutiveModels *> constitutiveModel = elements[0]->getConsModel();
                        dispBCManagement(nodes, DNod, ndim, ifs, npart, nodLocal, activeDofMapLocal,
                                         nActiveDof, ISactiveDof, duration, constitutiveModel, flagInstantaneous);
#ifdef PARALLEL
                        MPI_Barrier(MPI_COMM_WORLD); //Wait all the processors to reach this point
#endif
                    } else if (line2.compare(0, 33, "*BOUNDARY, TYPE=DISPLACEMENT INST") == 0) {
                        flagInstantaneous = true;

#ifdef PARALLEL
                        MPI_Barrier(MPI_COMM_WORLD); //Wait all the processors to reach this point
#endif
                        double duration = tf - t0;
                        vector<constitutiveModels *> constitutiveModel = elements[0]->getConsModel();
                        dispBCManagement(nodes, DNod, ndim, ifs, npart, nodLocal, activeDofMapLocal,
                                         nActiveDof, ISactiveDof, duration, constitutiveModel, flagInstantaneous);
#ifdef PARALLEL
                        MPI_Barrier(MPI_COMM_WORLD); //Wait all the processors to reach this point
#endif
                    } else if (line2.compare(0, 7, "*Forces") == 0) {
                        forceOutputs = new classPrintForces();
                        forceOutManagement(forceOutputs, ifs);
#ifdef PARALLEL
                        MPI_Barrier(MPI_COMM_WORLD); //Wait all the processors to reach this point
#endif
                    } else if (line2.compare(0, 8, "*Outputs") == 0) {
                        outputs = new classOutputs(); // Default outputs. U, VMS, Cauchy
                        outManagement(outputs, ifs);
#ifdef PARALLEL
                        MPI_Barrier(MPI_COMM_WORLD); //Wait all the processors to reach this point
#endif
                    } else if (line2.compare(0, 19, "*ADAPTIVE TIME STEP") == 0) {
                        getline(ifs, secondLine);
                        int tmp = int(secondLine.size());
                        int contc = 0;
                        commas.clear();
                        for (int j = 0; j < tmp; j++) {
                            if (secondLine.compare(j, 1, ",") == 0) {
                                commas.push_back(j);
                                contc++;
                            }
                        }
                        memset(tmpme2, 0, maxsize);
                        secondLine.copy(tmpme2, commas[0], 0);
                        
                        timeStepAdaptiveType.append(tmpme2);
                        for (int ii=1; ii< commas.size(); ii++)
                        {
                            memset(tmpme2, 0, maxsize);
                            secondLine.copy(tmpme2, commas[ii] - commas[ii-1] - 1, commas[ii-1] + 1);
                            timeStepAdaptiveData.push_back(atof(tmpme2));
                        }
                        memset(tmpme2, 0, maxsize);
                        secondLine.copy(tmpme2, secondLine.size()-commas[commas.size() - 1], commas[commas.size() - 1] + 1);
                        timeStepAdaptiveData.push_back(atof(tmpme2));
                    } else if (line2.compare(0, 5, "*Time") == 0) {
                        getline(ifs, secondLine);
                        int tmp = int(secondLine.size());
                        int contc = 0;
                        commas.clear();
                        for (int j = 0; j < tmp; j++) {
                            if (secondLine.compare(j, 1, ",") == 0) {
                                commas.push_back(j);
                                contc++;
                            }
                        }
                        memset(tmpme2, 0, maxsize);
                        secondLine.copy(tmpme2, commas[0], 0);
                        t0 = atof(tmpme2);
                        memset(tmpme2, 0, maxsize);
                        secondLine.copy(tmpme2, secondLine.size(), commas[1 - 1] + 1);
                        tf = atof(tmpme2);
                        secondLine.clear();
                        steppingPlan.setStartTime(t0);
                        steppingPlan.setEndTime(tf);
#ifdef PARALLEL
                        MPI_Barrier(MPI_COMM_WORLD); //Wait all the processors to reach this point
#endif
                    } else if (line2.compare(0, 19, "*TIME STEPPING PLAN") == 0) {
                        streampos oldposLoc = ifs.tellg();
                        while (getline(ifs, secondLine)) {
                            if (secondLine.compare(0, 1, "*") == 0) 
                            {
                                ifs.seekg(oldposLoc);
                                break;
                            } 
                            else
                            {
                                int tmp = int(secondLine.size());
                                int contc = 0;
                                commas.clear();
                                for (int j = 0; j < tmp; j++) {
                                    if (secondLine.compare(j, 1, ",") == 0) {
                                        commas.push_back(j);
                                        contc++;
                                    }
                                }
                                memset(tmpme2, 0, maxsize);
                                secondLine.copy(tmpme2, commas[0], 0);
                                double tt = atof(tmpme2);
                                memset(tmpme2, 0, maxsize);
                                secondLine.copy(tmpme2, secondLine.size(), commas[1 - 1] + 1);
                                int step = atof(tmpme2);
                                secondLine.clear();
                                steppingPlan.setNumberOfSteps(tt,step);
                                oldposLoc = ifs.tellg();
                            }
                        }
                        
                    } else if (line2.compare(0, 13, "*PETSC_SOLVER") == 0) {

                        /// Petsc solvers
                        std::vector<std::string> petsci_solvers{"richardson", "chebyshev", "cg", "groppcg", "pipecg",
                                                                "pipecgrr", "pipelcg", "pipeprcg", "pipecg2", "cgne",
                                                                "nash",
                                                                "stcg", "gltr", "nash", "stcg", "gltr", "fcg",
                                                                "pipefcg", "gmres", "pipefgmres", "fgmres", "lgmres",
                                                                "dgmres", "pgmres", "tcqmr",
                                                                "bcgs", "ibcgs", "fbcgs", "fbcgsr", "bcgsl", "pipebcgs",
                                                                "cgs", "tfqmr", "cr", "pipecr", "lsqr", "preonly",
                                                                "qcg", "bicg", "minres",
                                                                "symmlq", "lcd", "python", "gcr", "pipegcr", "tsirm",
                                                                "cgls", "fetidp", "hpddm", "lu", "cholesky"};


                        petsc_solver = line2.substr(15, 30);
                        if (std::find(petsci_solvers.begin(), petsci_solvers.end(), petsc_solver.c_str()) !=
                            petsci_solvers.end()) {
                            flag_solver = 1;
                        } else {
                            ERROR("PETSc solver is not properly defined");
                            exit(-1);
                        }

#ifdef PARALLEL
                        MPI_Barrier(MPI_COMM_WORLD); //Wait all the processors to reach this point
#endif

                    } else if (line2.compare(0, 9, "*PETSC_PC") == 0) {

                        /// Petsc preconditioners
                        std::vector<std::string> petsci_pc{"none", "jacobi", "sor", "shell", "bjacobi", "mg",
                                                           "eisenstat", "ilu", "icc", "asm", "gasm", "ksp", "composite",
                                                           "redundant", "spai",
                                                           "nn", "pbjacobi", "vpbjacobi", "mat", "hypre", "parms",
                                                           "fieldsplit", "tfs", "ml", "galerkin", "exotic", "cp",
                                                           "bfbt", "lsc",
                                                           "python", "pfmg", "syspfmg", "redistribute", "svd", "gamg",
                                                           "chowiluviennacl", "rowscalingviennacl", "saviennacl",
                                                           "bddc", "kaczmarz",
                                                           "telescope", "patch", "lmvm", "hmg", "deflation", "hpddm",
                                                           "hara"};

                        petsc_precon = line2.substr(11, 30);
                        if (std::find(petsci_pc.begin(), petsci_pc.end(), petsc_precon) != petsci_pc.end()) {
                            flag_pc = 1;
                        } else {
                            ERROR("PETSc preconditioner is not properly defined");
                            exit(-1);
                        }

#ifdef PARALLEL
                        MPI_Barrier(MPI_COMM_WORLD); //Wait all the processors to reach this point
#endif

                    }
                }
            }

#ifdef PARALLEL
            int rank;
            PetscErrorCode ierr;
            ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            if (rank == 0)
            {
                if (solverType != "EXPLICIT"){
                   if (flag_solver == 0) {
                        WARNING("PETSc solver is not defined by the user for the %s solver.  The LU solver is assigned by default", solverType.c_str());
                   } else if (flag_solver == 1) {
                        INFO("PETSc solver %s is assigned for the %s solver", petsc_solver.c_str(), solverType.c_str());
                   }
        
                   if (flag_pc == 0) {
                        WARNING("PETSc preconditioner is not defined by the user for the %s solver.", solverType.c_str());
                   } else if (flag_pc == 1) {
                       INFO("PETSc preconditioner %s is assigned for the %s solver", petsc_precon.c_str(), solverType.c_str());
                   }
                }
                if ((petsc_solver ==  "lu") || (petsc_solver == "cholesky")) // Direct solvers
                    WARNING("Please download the MUMPS solver, otherwise the direct solver will not work in parallel.");
            }
#else
            if (solverType != "EXPLICIT") {
                if (flag_solver == 0) {
                    WARNING("PETSc solver is not defined by the user for the %s solver. The LU solver is assigned by default",
                            solverType.c_str());
                } else if (flag_solver == 1) {
                    INFO("PETSc solver %s is assigned for the %s solver", petsc_solver.c_str(), solverType.c_str());
                }

                if (flag_pc == 0) {
                    WARNING("PETSc preconditioner is not defined by the user for the %s solver.", solverType.c_str());
                } else if (flag_pc == 1) {
                    INFO("PETSc preconditioner %s is assigned for the %s solver.", petsc_precon.c_str(),
                         solverType.c_str());
                }
            }
#endif

            if (solverType == "EXPLICIT") {
                sol = new classSolverExplicit(nodes, ndim, nDof, DNod, scaleFactor, outputs, t0, tf, solverType,
                                              activeDofMapLocal, nActiveDof, ISactiveDof, flag_AV, val_viscosity,
                                              forceOutputs, flagActivationSolver, petsc_solver, petsc_precon, flagFSI);
            } else if (solverType == "IMPLICIT") {
                sol = new classSolverImplicit(nodes, ndim, nDof, DNod, scaleFactor, outputs, t0, tf, solverType,
                                              activeDofMapLocal, nActiveDof, ISactiveDof, forceOutputs,
                                              flagActivationSolver, petsc_solver, petsc_precon, flagFSI);
            } else if (solverType == "IMPLICIT STATIC") {
                sol = new classSolverImplicitStatic(nodes, ndim, nDof, DNod, scaleFactor, outputs, t0, tf, solverType,
                                                    activeDofMapLocal, nActiveDof, ISactiveDof, forceOutputs,
                                                    flagActivationSolver, petsc_solver, petsc_precon, flagFSI);
            } else {
                ERROR("The solver %s is not supported by MuPhiSim", solverType.c_str());
                exit(-1);
            }
            if (numSteps >0)
            {
                sol->setNumberOfSteps(numSteps);
                numSteps = 0;
            };
            if (absTol >=0.)
            {
                sol->setAbsoluteTolerance(absTol);
                absTol = -1.;
            }
            if (relTol >=0.)
            {
                sol->setRelativeTolerance(relTol);
                relTol = -1;
            }
            if (maximalNumberIterations >0)
            {
                sol->setMaximalNumberOfIterations(maximalNumberIterations);
                maximalNumberIterations = 0;
            }
            if (stiffnessMatrixType != solvers::None)
            {
                sol->setStiffnessMatrixType(stiffnessMatrixType, fixedMatrixIteration);
                stiffnessMatrixType = solvers::None;
            }
            if (tangentByPerturbation)
            {
                sol->setTangentByPerturbation(true, perturbationTol);
                tangentByPerturbation = false;
                perturbationTol = 0;
            };
            
            // stepping plan
            sol->setTimeSteppingPlan(steppingPlan);
            steppingPlan.reset();
            
            // adaptive
            if (timeStepAdaptiveData.size() >0)
            {
                sol->setAdaptiveTimeStep(timeStepAdaptiveType.c_str(), timeStepAdaptiveData);
                timeStepAdaptiveType = "";
                timeStepAdaptiveData.clear();
            }
            
            if (maximalTimeStep >0)
            {
                sol->setMaximalTimeStep(maximalTimeStep);
                maximalTimeStep = -1.;
            }
            
            if (trueDispField.size() >0)
            {
                sol->allocateTrueDisplacementFieldForErrorAnalysis(trueDispField.c_str(),trueDispFieldData);
                trueDispField = "";
                trueDispFieldData.clear();
            };
            
            if (allowableNumFailedSteps >0)
            {
                sol->setAllowableFailedSteps(allowableNumFailedSteps);
                allowableNumFailedSteps=0;
            }
        
            //
            solver.push_back(sol);
            DNod.clear();
            activeDofMapLocal.clear();

            //nActiveDof=0;
            int nActiveDofSize = nActiveDof.size();
            for (int k = 0; k < nActiveDofSize; k++) {
                nActiveDof.pop_back();
            }

            ISactiveDof.clear();
        }
        solverType.clear();
        forceOutputs = nullptr;
        /// Reset the flags for the petsc solver and preconditioner
        flag_solver = 0;
        flag_pc = 0;
        petsc_precon = "none";
        petsc_solver = "lu";
    }

#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD); //Wait all the processors to reach this point
#endif
}

/*! \brief Load the force BCs from the input file
  @param[in] nodes Array of nodes
  @param[inout] NeumannBCs Array with the Neumann boundary conditions
  @param[in] ndim Dimension of the problem
  @param[in] ifs Pointer to the input file (.inp)
  @param[in] npart Processor ownership mapping of the nodes in the whole model npart[node i]=rank (empty for sequential simulations)
  @param[in] nodLocal Mapping of local nodes to global nodes numbering nodlocal[local node index]= global node index (empty for sequential simulations)
*/
extern void forceBCManagement(vector<classNodes *> &nodes, vector<classNeumannBCs *> &NeumannBCs, int ndim,
                              ifstream &ifs, vector<int> &npart, vector<int> &nodLocal) {

    vector<classNeumann *> NeumannNodes;
    string type = "FORCE";
    vector<double> tempxyz;
    vector<int> commas;
    classNeumann *tempNeumann;
    classNeumannBCs *nodForces;
    int maxsize = 100;
    char tmpme2[maxsize], tempU[maxsize];
    memset(tmpme2, 0, maxsize);
    memset(tempU, 0, maxsize);
    double tmp = 0, columns = 3;
    string line;
    int tempme = 0;
    streampos oldpos;

#ifdef PARALLEL
    vector<int>::iterator it;
    int rank;
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    // Time distribution
    getline(ifs, line);
    int contc = 0;
    double t0 = 0, tf = 0, Ff = 0;
    for (int j = 0; j < line.size(); j++) {
        if (line.compare(j, 1, ",") == 0) {
            commas.push_back(j);
            contc++;
        }
    }
    line.copy(tmpme2, commas[0], 0);
    t0 = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[1] - commas[1 - 1] - 1, commas[1 - 1] + 1);
    tf = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, line.size(), commas[2 - 1] + 1);
    Ff = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.clear();
    commas.clear();

    while (getline(ifs, line)) {
        if (line.compare(0, 1, "*") == 0) {
            ifs.seekg(oldpos);
            break;
        } else {
            tmp = int(line.size());
            contc = 0;
            vector<double> F(columns, 0);
            vector<double> Fff(columns, 0);
            for (int j = 0; j < tmp; j++) {
                if (line.compare(j, 1, ",") == 0) {
                    commas.push_back(j);
                    contc++;
                }
            }
            line.copy(tmpme2, commas[0], 0);
            for (int k = 1; k < columns; k++) {
                line.copy(tempU, commas[k] - commas[k - 1] - 1, commas[k - 1] + 1);
                F[k - 1] = atof(tempU);
                memset(tempU, 0, maxsize);
            }
            line.copy(tempU, tmp - commas[columns - 1] - 1, commas[columns - 1] + 1);
            F[2] = atof(tempU);
            tempme = atoi(tmpme2) - 1;
            for (int i = 0; i < ndim; i++) {
                Fff[i] = F[i] * Ff;
            }

#ifdef PARALLEL
            if (npart[tempme] == rank) {
                it = find(nodLocal.begin(), nodLocal.end(), tempme);
                tempxyz = nodes[(it - nodLocal.begin())]->getXYZ();
                tempNeumann = new classNeumann(ndim, (it - nodLocal.begin()), tempxyz, Fff, "Linear", t0, tf);
                NeumannNodes.push_back(tempNeumann);
                //it.clear();
            }
#else
            tempxyz = nodes[tempme]->getXYZ();
            tempNeumann = new classNeumann(ndim, tempme, tempxyz, Fff, "Linear", t0, tf);
            NeumannNodes.push_back(tempNeumann);
#endif


            commas.clear();
            F[0] = 0;
            F[1] = 0;
            F[2] = 0;
            memset(tmpme2, 0, maxsize);
            memset(tempU, 0, maxsize);
        }
        oldpos = ifs.tellg();
    }
    nodForces = new classNodalForces(type, NeumannNodes);
    NeumannBCs.push_back(nodForces);

}


/*! \brief Load the disp BCs from the input file
  @param[in] nodes Array of nodes
  @param[inout] DirichletNodes Array with all dirichlet nodes
  @param[in] ndim Dimension of the problem
  @param[in] ifs Pointer to the input file (.inp)
  @param[in] npart Processor ownership mapping of the nodes in the whole model npart[node i]=rank (empty for sequential simulations)
  @param[out] nodLocal Mapping of local nodes to global nodes numbering nodlocal[local node index]= global node index (empty for sequential simulations)
  @param[in] activeDofMapLocal Array for parallel simulations with the local active Dof of a processor (empty in sequential simulations)
  @param[in] nActiveDof Global number of active Dof in parallel simulations (zero in sequential simulations)
  @param[in] ISactiveDof Array for parallel simulations with the mapping for the global displacements vector (empty in sequential simulations)
  @param[in] duration double for recording the time duration of BC

*/
//displacement management to take into account of input time for calculating correct loading velocity, and correcting the global to active mapping
extern void
dispBCManagement(vector<classNodes *> &nodes, vector<classDirichlet *> &DirichletNodes, int ndim, ifstream &ifs,
                 vector<int> &npart, vector<int> &nodLocal, vector<int> &activeDofMapLocal, vector<int> &nActiveDof,
                 vector<int> &ISactiveDof, double duration, vector<constitutiveModels *> constitutiveModel,
                 bool flagInstantaneous) {

    string type = "DISPLACEMENT";
    classDirichlet *tempDirichlet;
    int numNodes = nodes.size();
    int maxsize = 100;
    char tmpme2[maxsize], tempU[maxsize], tempDirection[maxsize];
    memset(tmpme2, 0, maxsize);
    memset(tempU, 0, maxsize);
    memset(tempDirection, 0, maxsize);
    string line;
    int tmp = 0;
    double UU = 0;
    vector<double> tempxyz;
    vector<int> commas;
    int tempme = 0;
    streampos oldpos;

    int nbrDofExtra = 0;
    for (int i = 1; i < constitutiveModel.size(); i++) { // i start from 1 in order to don't take into account the mechanical model
        nbrDofExtra += constitutiveModel[i]->getNbrDofConsMod();
    }
    int nbrDofStochastic = constitutiveModel[0]->getNbrDofConsMod() - ndim;
    
    vector<int>::iterator it;

#ifdef PARALLEL
    vector<int> auxActDof((ndim + nbrDofExtra + nbrDofStochastic) * npart.size(), 0);
    vector<int> activeDofMap((ndim + nbrDofExtra + nbrDofStochastic) * npart.size(), 0);
    vector<int> DirichletNodesGlobal;
    vector<int> DirichletNodesGlobalExtra;
    activeDofMapLocal.resize((ndim + nbrDofExtra + nbrDofStochastic) * nodLocal.size(), 0);
    int numGlobalNodes = npart.size();
    PetscErrorCode ierr;
    int rank = 0;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    CHKERRV(ierr);
    MPI_Barrier(MPI_COMM_WORLD);
    //#endif
#else
    vector<int> auxActDof((ndim + nbrDofExtra + nbrDofStochastic) * numNodes, 0);
    vector<int> activeDofMap((ndim + nbrDofExtra + nbrDofStochastic) * numNodes, 0);
    vector<int> DirichletNodesGlobal;
    vector<int> DirichletNodesGlobalExtra;
#endif

#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD); //Wait all the processors to reach this point
#endif
    while (getline(ifs, line)) {
        if (line.compare(0, 1, "*") == 0) {
            ifs.seekg(oldpos);
            break;
        } else {
            tmp = int(line.size());
            int contc = 0;
            commas.clear();
            for (int j = 0; j < line.size(); j++) {
                if (line.compare(j, 1, ",") == 0) {
                    commas.push_back(j);
                    contc++;
                }
            }
            line.copy(tmpme2, commas[0], 0);
            // What direction
            // 1-x, 2-y, 3-z
            line.copy(tempDirection, commas[1] - commas[1 - 1] - 1, commas[1 - 1] + 1);
            line.copy(tempU, tmp - commas[2 - 1] - 1, commas[2 - 1] + 1);
            tempme = atoi(tmpme2) - 1;
            if (npart.size() == 0) {
                int direction = atoi(tempDirection);
                UU = atof(tempU);
                tempxyz = nodes[tempme]->getXYZ();

                /*if(direction>=3){ Kamel
                    ERROR("The direction %i at node %i does not make any sense", direction, tempme+1);
                    INFO("Reminder: 0-x, 1-y and 2-z");
                    exit(-1);
                }*/

                if (direction >= (ndim + nbrDofExtra + nbrDofStochastic)) { //Kamel
                    ERROR("The direction %i at node %i has not been recognized", direction, tempme + 1);
                    exit(-1);
                }
                if (flagInstantaneous) {
                    //if flagInstantaneous = true Dirichlet BC are applied instantaneously
                    tempDirichlet = new classDirichlet(ndim, tempme, tempxyz, UU, direction, numNodes);
                } else {
                    //if flagInstantaneous = false Dirichlet BC they are applied as a ramp
                    tempDirichlet = new classDirichlet(ndim, tempme, tempxyz, UU, direction, numNodes, duration);
                }
                DirichletNodes.push_back(tempDirichlet);
                //test code
                if (direction < ndim) {
                    DirichletNodesGlobal.push_back(((direction) * numNodes) + tempme);
                } else {
                    DirichletNodesGlobalExtra.push_back(((direction) * numNodes) + tempme);
                }
            } else { //parallel

#ifdef PARALLEL
                int direction = atoi(tempDirection);
                if (direction >= (ndim + nbrDofExtra + nbrDofStochastic)) { //Kamel
                    ERROR("The direction %i at node %i has not been recognized", direction, tempme + 1);
                    exit(-1);
                }
                if (direction < ndim) {
                    DirichletNodesGlobal.push_back(((direction) * numGlobalNodes) + tempme);
                } else {
                    DirichletNodesGlobalExtra.push_back(((direction) * numGlobalNodes) + tempme);
                }

                if (npart[tempme] == rank) {
                    it = find(nodLocal.begin(), nodLocal.end(), tempme);
                    tempme = it - nodLocal.begin();
                    UU = atof(tempU);
                    tempxyz = nodes[tempme]->getXYZ();

                    /*if(direction>=3){
                        ERROR("The direction %i at node %i does not make any sense", direction, tempme+1);
                        INFO("Reminder: 0-x, 1-y and 2-z");
                        exit(-1);
                    }*/

                    if (flagInstantaneous) {
                        //if flagInstantaneous = true Dirichlet BC are applied instantaneously
                        tempDirichlet = new classDirichlet(ndim, tempme, tempxyz, UU, direction, numNodes);
                    } else {
                        //if flagInstantaneous = false Dirichlet BC they are applied as a ramp
                        tempDirichlet = new classDirichlet(ndim, tempme, tempxyz, UU, direction, numNodes, duration);
                    }

                    DirichletNodes.push_back(tempDirichlet);

                }
#endif
            }
            memset(tmpme2, 0, maxsize);
            memset(tempDirection, 0, maxsize);
            memset(tempU, 0, maxsize);
            commas.clear();
        }
        oldpos = ifs.tellg();
    }

#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD); //Wait all the processors to reach this point
#endif

#ifdef PARALLEL
    if (!(npart.size() == 0)) {

        sort(DirichletNodesGlobal.begin(), DirichletNodesGlobal.end());
        nActiveDof.push_back(npart.size() * ndim - DirichletNodesGlobal.size());

        sort(DirichletNodesGlobalExtra.begin(), DirichletNodesGlobalExtra.end());
        nActiveDof.push_back(npart.size() * (nbrDofExtra + nbrDofStochastic) - DirichletNodesGlobalExtra.size());

        //Building up the mapping vector
        //First filling the substraction dirichlet indexes
        int mark = 0;
        /*Kamel
        for(int k=0;k<DirichletNodesGlobal.size();k++){
            int mark2=DirichletNodesGlobal[k];
            fill(auxActDof.begin()+(mark),auxActDof.begin()+mark2,-k);
            mark=DirichletNodesGlobal[k];
        }
        fill(auxActDof.begin()+DirichletNodesGlobal[DirichletNodesGlobal.size()-1],auxActDof.end(),-(DirichletNodesGlobal.size()));
        */
        // ****** Kamel *******
        int mark2 = 0;
        for (int k = 0; k < (DirichletNodesGlobal.size() + DirichletNodesGlobalExtra.size()); k++) {
            if (k < DirichletNodesGlobal.size()) {
                mark2 = DirichletNodesGlobal[k];
            } else {
                mark2 = DirichletNodesGlobalExtra[k - DirichletNodesGlobal.size()];
            }
            fill(auxActDof.begin() + (mark), auxActDof.begin() + mark2, -k);

            if (k < DirichletNodesGlobal.size()) {
                mark = DirichletNodesGlobal[k];
            } else {
                mark = DirichletNodesGlobalExtra[k - DirichletNodesGlobal.size()];
            }
        }
        if (DirichletNodesGlobalExtra.size() == 0) {
            fill(auxActDof.begin() + DirichletNodesGlobal[DirichletNodesGlobal.size() - 1], auxActDof.end(),
                 -(DirichletNodesGlobal.size()));
        } else {
            fill(auxActDof.begin() + DirichletNodesGlobalExtra[DirichletNodesGlobalExtra.size() - 1], auxActDof.end(),
                 -(DirichletNodesGlobal.size() + DirichletNodesGlobalExtra.size()));
        }
        // ****** Kamel *******

        //Now creating the global activeDofMap
        for (int j = 0; j < activeDofMap.size(); j++) {
            activeDofMap[j] = j + auxActDof[j];
        }

        //Forcing Dirichlet DOF to be -1
        for (int j = 0; j < DirichletNodesGlobal.size(); j++) {
            activeDofMap[DirichletNodesGlobal[j]] = -1;
        }
        for (int j = 0; j < DirichletNodesGlobalExtra.size(); j++) {
            activeDofMap[DirichletNodesGlobalExtra[j]] = -1;
        }
        for (int j = 0; j < activeDofMap.size(); j++) {
            if (!(activeDofMap[j] == -1)) {
                ISactiveDof.push_back(j);
            }
        }
        //Recovering the local nodes with the global active numbering
        for (int j = 0; j < nodLocal.size(); j++) {
            for (int i = 0; i < (ndim + nbrDofExtra + nbrDofStochastic); i++) {
                activeDofMapLocal[(nodLocal.size() * i) + j] = activeDofMap[(i * npart.size()) + nodLocal[j]];
            }
        }
        activeDofMap.clear();
        DirichletNodesGlobal.clear();
    }
#else
    //store the active Dof info to class varibles in sequential mode as well
    sort(DirichletNodesGlobal.begin(), DirichletNodesGlobal.end());
    nActiveDof.push_back(numNodes * ndim - DirichletNodesGlobal.size());

    sort(DirichletNodesGlobalExtra.begin(), DirichletNodesGlobalExtra.end());
    nActiveDof.push_back(numNodes * (nbrDofExtra + nbrDofStochastic) - DirichletNodesGlobalExtra.size());

    int count = 0;
    for (int i = 0; i < activeDofMap.size(); i++) {
        int nextDirichlet = -1; // Just to be sure that the next test "i == nextDirichlet" remains false if DirichletNodesGlobalExtra.size() == 0
        if (i <= DirichletNodesGlobal[DirichletNodesGlobal.size() - 1]) {
            nextDirichlet = DirichletNodesGlobal[count];
        } else {
            if (DirichletNodesGlobalExtra.size() != 0) {
                if (i <= DirichletNodesGlobalExtra[DirichletNodesGlobalExtra.size() -
                                                   1]) { //without this line, the code is out of size
                    nextDirichlet = DirichletNodesGlobalExtra[count - DirichletNodesGlobal.size()];
                }
            }
        }
        if (i == nextDirichlet) {
            activeDofMap[i] = -1;
            count++;
        } else {
            activeDofMap[i] = i - count;
            ISactiveDof.push_back(i);
        }
    }
    activeDofMapLocal = activeDofMap;
    activeDofMap.clear();
    DirichletNodesGlobal.clear();
#endif

}

/*! \brief Load the force BCs from the input file
  @param[in] elements Array with all elements (or background cells) in the domain
  @param[in] nodes Array of nodes
  @param[inout] NeumannBCs Array with the Neumann boundary conditions
  @param[in] ndim Dimension of the problem
  @param[in] ifs Pointer to the input file (.inp)
  @param[in] npart Processor ownership mapping of the nodes in the whole model npart[node i]=rank (empty for sequential simulations)
  @param[in] nodLocal Mapping of local nodes to global nodes numbering nodlocal[local node index]= global node index (empty for sequential simulations)
*/
extern void tractionInstBCManagement(vector<classElements *> &elements, vector<classNodes *> &nodes,
                                     vector<classNeumannBCs *> &NeumannBCs, int ndim, ifstream &ifs, vector<int> &epart,
                                     vector<int> &elemLocal, int surface_order) {

    classNeumannBCs *nodForces;
    string type = "TRACTION INST";
    vector<double> tempxyz;
    vector<int> commas;
    int maxsize = 100;
    char tmpme2[maxsize], tempU[maxsize];
    memset(tmpme2, 0, maxsize);
    memset(tempU, 0, maxsize);
    double tmp = 0;
    string line;
    double Tx = 0;
    double Ty = 0;
    double Tz = 0;
    int S = 0;
    vector<classElements *> elemsPress;
    streampos oldpos;

    getline(ifs, line);
    tmp = int(line.size());
    int contc = 0;
    for (int j = 0; j < tmp; j++) {
        if (line.compare(j, 1, ",") == 0) {
            commas.push_back(j);
            contc++;
        }
    }

    double t0 = 0;
    double tf = 0;
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[0], 0);
    t0 = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[1] - commas[1 - 1] - 1, commas[1 - 1] + 1);
    tf = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[2] - commas[2 - 1] - 1, commas[2 - 1] + 1);
    Tx = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[3] - commas[3 - 1] - 1, commas[3 - 1] + 1);
    Ty = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[4] - commas[4 - 1] - 1, commas[4 - 1] + 1);
    Tz = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, line.size(), commas[5 - 1] + 1);
    S = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.clear();
#ifdef PARALLEL
    vector<int>::iterator it;
    int rank;
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    while (getline(ifs, line)) {
        if (line.compare(0, 1, "*") == 0) {
            ifs.seekg(oldpos);
            break;
        } else {
            //PARALLEL
            if (epart.size() != 0) {
#ifdef PARALLEL
                int size = int(line.size());
                line.copy(tmpme2, size, 0);
                vector<classElements *> subElements;
                int mainEle = atoi(tmpme2) - 1;
                memset(tmpme2, 0, maxsize);
                if (epart[mainEle] == rank) {
                    it = find(elemLocal.begin(), elemLocal.end(), mainEle);
                    mainEle = it - elemLocal.begin();
                    elements[mainEle]->getSubElement(elemsPress.size(), S, subElements, nodes, mainEle,
                                                     surface_order); //is me a unique identifier? dont think it is being used or matters
                    for (int i = 0; i < subElements.size(); i++) {
                        elemsPress.push_back(subElements[i]);
                    }

                }
#endif
                //SERIAL
            } else {
                int size = int(line.size());
                line.copy(tmpme2, size, 0);
                vector<classElements *> subElements;
                int mainEle = atoi(tmpme2) - 1;
                memset(tmpme2, 0, maxsize);
                elements[mainEle]->getSubElement(elemsPress.size(), S, subElements, nodes, mainEle,
                                                 surface_order); //is me a unique identifier? dont think it is being used or matters
                for (int i = 0; i < subElements.size(); i++) {
                    elemsPress.push_back(subElements[i]);
                }
            }
        }
        oldpos = ifs.tellg();
    }
    nodForces = new classTractionInst(type, elemsPress, Tx, Ty, Tz, t0, tf);
    NeumannBCs.push_back(nodForces);
}

/*! \brief Load the force BCs from the input file
  @param[in] elements Array with all elements (or background cells) in the domain
  @param[in] nodes Array of nodes
  @param[inout] NeumannBCs Array with the Neumann boundary conditions
  @param[in] ndim Dimension of the problem
  @param[in] ifs Pointer to the input file (.inp)
  @param[in] npart Processor ownership mapping of the nodes in the whole model npart[node i]=rank (empty for sequential simulations)
  @param[in] nodLocal Mapping of local nodes to global nodes numbering nodlocal[local node index]= global node index (empty for sequential simulations)
*/
extern void tractionRampBCManagement(vector<classElements *> &elements, vector<classNodes *> &nodes,
                                     vector<classNeumannBCs *> &NeumannBCs, int ndim, ifstream &ifs, vector<int> &epart,
                                     vector<int> &elemLocal, int surface_order) {

    classNeumannBCs *nodForces;
    string type = "TRACTION RAMP";
    vector<double> tempxyz;
    vector<int> commas;
    int maxsize = 100;
    char tmpme2[maxsize], tempU[maxsize];
    memset(tmpme2, 0, maxsize);
    memset(tempU, 0, maxsize);
    double tmp = 0;
    string line;
    double Tx = 0;
    double Ty = 0;
    double Tz = 0;
    int S = 0;
    vector<classElements *> elemsPress;
    streampos oldpos;

    getline(ifs, line);
    tmp = int(line.size());
    int contc = 0;
    for (int j = 0; j < tmp; j++) {
        if (line.compare(j, 1, ",") == 0) {
            commas.push_back(j);
            contc++;
        }
    }

    double t0 = 0;
    double tf = 0;
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[0], 0);
    t0 = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[1] - commas[1 - 1] - 1, commas[1 - 1] + 1);
    tf = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[2] - commas[2 - 1] - 1, commas[2 - 1] + 1);
    Tx = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[3] - commas[3 - 1] - 1, commas[3 - 1] + 1);
    Ty = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[4] - commas[4 - 1] - 1, commas[4 - 1] + 1);
    Tz = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, line.size(), commas[5 - 1] + 1);
    S = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.clear();
#ifdef PARALLEL
    vector<int>::iterator it;
    int rank;
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    while (getline(ifs, line)) {
        if (line.compare(0, 1, "*") == 0) {
            ifs.seekg(oldpos);
            break;
        } else {
            //PARALLEL
            if (epart.size() != 0) {
#ifdef PARALLEL
                int size = int(line.size());
                line.copy(tmpme2, size, 0);
                vector<classElements *> subElements;
                int mainEle = atoi(tmpme2) - 1;
                memset(tmpme2, 0, maxsize);
                if (epart[mainEle] == rank) {
                    it = find(elemLocal.begin(), elemLocal.end(), mainEle);
                    mainEle = it - elemLocal.begin();
                    elements[mainEle]->getSubElement(elemsPress.size(), S, subElements, nodes, mainEle,
                                                     surface_order); //is me a unique identifier? dont think it is being used or matters
                    for (int i = 0; i < subElements.size(); i++) {
                        elemsPress.push_back(subElements[i]);
                    }

                }
#endif
                //SERIAL
            } else {
                int size = int(line.size());
                line.copy(tmpme2, size, 0);
                vector<classElements *> subElements;
                int mainEle = atoi(tmpme2) - 1;
                memset(tmpme2, 0, maxsize);
                elements[mainEle]->getSubElement(elemsPress.size(), S, subElements, nodes, mainEle,
                                                 surface_order); //is me a unique identifier? dont think it is being used or matters
                for (int i = 0; i < subElements.size(); i++) {
                    elemsPress.push_back(subElements[i]);
                }
            }
        }
        oldpos = ifs.tellg();
    }
    nodForces = new classTractionRamp(type, elemsPress, Tx, Ty, Tz, t0, tf);
    NeumannBCs.push_back(nodForces);
}

/*! \brief Load the stochastic force BCs from the input file
  @param[in] elements Array with all elements (or background cells) in the domain
  @param[in] nodes Array of nodes
  @param[inout] NeumannBCs Array with the Neumann boundary conditions
  @param[in] ndim Dimension of the problem
  @param[in] ifs Pointer to the input file (.inp)
  @param[in] npart Processor ownership mapping of the nodes in the whole model npart[node i]=rank (empty for sequential simulations)
  @param[in] nodLocal Mapping of local nodes to global nodes numbering nodlocal[local node index]= global node index (empty for sequential simulations)
*/
extern void pressureStochasticRampBCManagement(vector<classElements *> &elements, vector<classNodes *> &nodes,
                                     vector<classNeumannBCs *> &NeumannBCs, int ndim, ifstream &ifs, vector<int> &epart,
                                     vector<int> &elemLocal, int surface_order) {

    classNeumannBCs *nodForces;
    string type = "STOCHASTIC PRESSURE RAMP";
    vector<double> tempxyz;
    vector<int> commas;
    int maxsize = 100;
    char tmpme2[maxsize], tempU[maxsize];
    memset(tmpme2, 0, maxsize);
    memset(tempU, 0, maxsize);
    double tmp = 0;
    string line;
    double P = 0;
    double P_var = 0;
    int S = 0;
    vector<classElements *> elemsPress;
    streampos oldpos;

    getline(ifs, line);
    tmp = int(line.size());
    int contc = 0;
    for (int j = 0; j < tmp; j++) {
        if (line.compare(j, 1, ",") == 0) {
            commas.push_back(j);
            contc++;
        }
    }

    double t0 = 0;
    double tf = 0;
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[0], 0);
    t0 = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[1] - commas[1 - 1] - 1, commas[1 - 1] + 1);
    tf = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[2] - commas[2 - 1] - 1, commas[2 - 1] + 1);
    P = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[3] - commas[3 - 1] - 1, commas[3 - 1] + 1);
    P_var = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, line.size(), commas[4 - 1] + 1);
    S = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.clear();
#ifdef PARALLEL
    vector<int>::iterator it;
    int rank;
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    while (getline(ifs, line)) {
        if (line.compare(0, 1, "*") == 0) {
            ifs.seekg(oldpos);
            break;
        } else {
            //PARALLEL
            if (epart.size() != 0) {
#ifdef PARALLEL
                int size = int(line.size());
                line.copy(tmpme2, size, 0);
                vector<classElements *> subElements;
                int mainEle = atoi(tmpme2) - 1;
                memset(tmpme2, 0, maxsize);
                if (epart[mainEle] == rank) {
                    it = find(elemLocal.begin(), elemLocal.end(), mainEle);
                    mainEle = it - elemLocal.begin();
                    elements[mainEle]->getSubElement(elemsPress.size(), S, subElements, nodes, mainEle,
                                                     surface_order); //is me a unique identifier? dont think it is being used or matters
                    for (int i = 0; i < subElements.size(); i++) {
                        elemsPress.push_back(subElements[i]);
                    }

                }
#endif
                //SERIAL
            } else {
                int size = int(line.size());
                line.copy(tmpme2, size, 0);
                vector<classElements *> subElements;
                int mainEle = atoi(tmpme2) - 1;
                memset(tmpme2, 0, maxsize);
                elements[mainEle]->getSubElement(elemsPress.size(), S, subElements, nodes, mainEle,
                                                 surface_order); //is me a unique identifier? dont think it is being used or matters
                for (int i = 0; i < subElements.size(); i++) {
                    elemsPress.push_back(subElements[i]);
                }
            }
        }
        oldpos = ifs.tellg();
    }

    vector<constitutiveModels *> elementsCons = elements[0]->getConsModel();
    int nbr_dof;
    nbr_dof = (elementsCons[0]->getNbrDofConsMod());
    int nstoch;
    nstoch = nbr_dof/ndim;
    vector<double> C(nstoch*nstoch*nstoch);
        const classStochasticHyperElasticStVenantKirchhoff* firstLaw = dynamic_cast<const classStochasticHyperElasticStVenantKirchhoff*>(elementsCons[0]);
    if (firstLaw == NULL)
    {
        ERROR("classStochasticHyperElasticStVenantKirchhoff must be used");
        exit(-1);
    }
    C = firstLaw->getC_Stoch();

    nodForces = new classStochasticPressureRamp(type, elemsPress, P, P_var, t0, tf, nstoch, C);
    NeumannBCs.push_back(nodForces);

}

/*! \brief Load the force BCs from the input file
  @param[in] elements Array with all elements (or background cells) in the domain
  @param[in] nodes Array of nodes
  @param[inout] NeumannBCs Array with the Neumann boundary conditions
  @param[in] ndim Dimension of the problem
  @param[in] ifs Pointer to the input file (.inp)
  @param[in] npart Processor ownership mapping of the nodes in the whole model npart[node i]=rank (empty for sequential simulations)
  @param[in] nodLocal Mapping of local nodes to global nodes numbering nodlocal[local node index]= global node index (empty for sequential simulations)
*/
extern void pressureRampBCManagement(vector<classElements *> &elements, vector<classNodes *> &nodes,
                                     vector<classNeumannBCs *> &NeumannBCs, int ndim, ifstream &ifs, vector<int> &epart,
                                     vector<int> &elemLocal, int surface_order) {

    classNeumannBCs *nodForces;
    string type = "PRESSURE RAMP";
    vector<double> tempxyz;
    vector<int> commas;
    int maxsize = 100;
    char tmpme2[maxsize], tempU[maxsize];
    memset(tmpme2, 0, maxsize);
    memset(tempU, 0, maxsize);
    double tmp = 0;
    string line;
    double P = 0;
    int S = 0;
    vector<classElements *> elemsPress;
    streampos oldpos;

    getline(ifs, line);
    tmp = int(line.size());
    int contc = 0;
    for (int j = 0; j < tmp; j++) {
        if (line.compare(j, 1, ",") == 0) {
            commas.push_back(j);
            contc++;
        }
    }

    double t0 = 0;
    double tf = 0;
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[0], 0);
    t0 = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[1] - commas[1 - 1] - 1, commas[1 - 1] + 1);
    tf = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[2] - commas[2 - 1] - 1, commas[2 - 1] + 1);
    P = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, line.size(), commas[3 - 1] + 1);
    S = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.clear();
#ifdef PARALLEL
    vector<int>::iterator it;
    int rank;
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    while (getline(ifs, line)) {
        if (line.compare(0, 1, "*") == 0) {
            ifs.seekg(oldpos);
            break;
        } else {
            //PARALLEL
            if (epart.size() != 0) {
#ifdef PARALLEL
                int size = int(line.size());
                line.copy(tmpme2, size, 0);
                vector<classElements *> subElements;
                int mainEle = atoi(tmpme2) - 1;
                memset(tmpme2, 0, maxsize);
                if (epart[mainEle] == rank) {
                    it = find(elemLocal.begin(), elemLocal.end(), mainEle);
                    mainEle = it - elemLocal.begin();
                    elements[mainEle]->getSubElement(elemsPress.size(), S, subElements, nodes, mainEle,
                                                     surface_order); //is me a unique identifier? dont think it is being used or matters
                    for (int i = 0; i < subElements.size(); i++) {
                        elemsPress.push_back(subElements[i]);
                    }

                }
#endif
                //SERIAL
            } else {
                int size = int(line.size());
                line.copy(tmpme2, size, 0);
                vector<classElements *> subElements;
                int mainEle = atoi(tmpme2) - 1;
                memset(tmpme2, 0, maxsize);
                elements[mainEle]->getSubElement(elemsPress.size(), S, subElements, nodes, mainEle,
                                                 surface_order); //is me a unique identifier? dont think it is being used or matters
                for (int i = 0; i < subElements.size(); i++) {
                    elemsPress.push_back(subElements[i]);
                }
            }
        }
        oldpos = ifs.tellg();
    }
    nodForces = new classPressureRamp(type, elemsPress, P, t0, tf);
    NeumannBCs.push_back(nodForces);

}

/*! \brief Load the force BCs from the input file
  @param[in] elements Array with all elements (or background cells) in the domain
  @param[in] nodes Array of nodes
  @param[inout] NeumannBCs Array with the Neumann boundary conditions
  @param[in] ndim Dimension of the problem
  @param[in] ifs Pointer to the input file (.inp)
  @param[in] npart Processor ownership mapping of the nodes in the whole model npart[node i]=rank (empty for sequential simulations)
  @param[in] nodLocal Mapping of local nodes to global nodes numbering nodlocal[local node index]= global node index (empty for sequential simulations)
*/
extern void pressureInstBCManagement(vector<classElements *> &elements, vector<classNodes *> &nodes,
                                     vector<classNeumannBCs *> &NeumannBCs, int ndim, ifstream &ifs, vector<int> &epart,
                                     vector<int> &elemLocal, int surface_order) {

    classNeumannBCs *nodForces;
    string type = "PRESSURE INST";
    vector<double> tempxyz;
    vector<int> commas;
    int maxsize = 100;
    char tmpme2[maxsize], tempU[maxsize];
    memset(tmpme2, 0, maxsize);
    memset(tempU, 0, maxsize);
    double tmp = 0;
    string line;
    double P = 0;
    int S = 0;
    vector<classElements *> elemsPress;
    streampos oldpos;

    getline(ifs, line);
    tmp = int(line.size());
    int contc = 0;
    for (int j = 0; j < tmp; j++) {
        if (line.compare(j, 1, ",") == 0) {
            commas.push_back(j);
            contc++;
        }
    }

    double t0 = 0;
    double tf = 0;
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[0], 0);
    t0 = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[1] - commas[1 - 1] - 1, commas[1 - 1] + 1);
    tf = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[2] - commas[2 - 1] - 1, commas[2 - 1] + 1);
    P = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, line.size(), commas[3 - 1] + 1);
    S = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.clear();
#ifdef PARALLEL
    vector<int>::iterator it;
    int rank;
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    while (getline(ifs, line)) {
        if (line.compare(0, 1, "*") == 0) {
            ifs.seekg(oldpos);
            break;
        } else {
            //PARALLEL
            if (epart.size() != 0) {
#ifdef PARALLEL
                int size = int(line.size());
                line.copy(tmpme2, size, 0);
                vector<classElements *> subElements;
                int mainEle = atoi(tmpme2) - 1;
                memset(tmpme2, 0, maxsize);
                if (epart[mainEle] == rank) {
                    it = find(elemLocal.begin(), elemLocal.end(), mainEle);
                    mainEle = it - elemLocal.begin();
                    elements[mainEle]->getSubElement(elemsPress.size(), S, subElements, nodes, mainEle,
                                                     surface_order); //is me a unique identifier? dont think it is being used or matters
                    for (int i = 0; i < subElements.size(); i++) {
                        elemsPress.push_back(subElements[i]);
                    }
                }
#endif
                //SERIAL
            } else {
                int size = int(line.size());
                line.copy(tmpme2, size, 0);
                vector<classElements *> subElements;
                int mainEle = atoi(tmpme2) - 1;
                memset(tmpme2, 0, maxsize);
                elements[mainEle]->getSubElement(elemsPress.size(), S, subElements, nodes, mainEle,
                                                 surface_order); //is me a unique identifier? dont think it is being used or matters
                for (int i = 0; i < subElements.size(); i++) {
                    elemsPress.push_back(subElements[i]);
                }
            }
        }
        oldpos = ifs.tellg();
    }
    nodForces = new classPressureInst(type, elemsPress, P, t0, tf);
    NeumannBCs.push_back(nodForces);

}

/*! \brief Load the force BCs from the input file
  @param[in] elements Array with all elements (or background cells) in the domain
  @param[in] nodes Array of nodes
  @param[inout] NeumannBCs Array with the Neumann boundary conditions
  @param[in] ndim Dimension of the problem
  @param[in] ifs Pointer to the input file (.inp)
  @param[in] npart Processor ownership mapping of the nodes in the whole model npart[node i]=rank (empty for sequential simulations)
  @param[in] nodLocal Mapping of local nodes to global nodes numbering nodlocal[local node index]= global node index (empty for sequential simulations)
*/

extern void hertzianRampBCManagement(vector<classElements *> &elements, vector<classNodes *> &nodes,
                                     vector<classNeumannBCs *> &NeumannBCs, int ndim, ifstream &ifs, vector<int> &epart,
                                     vector<int> &elemLocal, int surface_order) {

    classNeumannBCs *nodForces;
    string type = "HERTZIAN RAMP";
    vector<double> tempxyz;
    vector<int> commas;
    int maxsize = 100;
    char tmpme2[maxsize], tempU[maxsize];
    memset(tmpme2, 0, maxsize);
    memset(tempU, 0, maxsize);
    double tmp = 0;
    string line;
    double P = 0;
    double T = 0;
    double mu;
    double a;
    int S = 0;
    vector<double> center;
    vector<classElements *> elemsPress;
    streampos oldpos;

    getline(ifs, line);
    tmp = int(line.size());
    int contc = 0;
    for (int j = 0; j < tmp; j++) {
        if (line.compare(j, 1, ",") == 0) {
            commas.push_back(j);
            contc++;
        }
    }

    double t0 = 0;
    double tf = 0;
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[0], 0);
    t0 = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[1] - commas[1 - 1] - 1, commas[1 - 1] + 1);
    tf = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[2] - commas[2 - 1] - 1, commas[2 - 1] + 1);
    P = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[3] - commas[3 - 1] - 1, commas[3 - 1] + 1);
    T = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[4] - commas[4 - 1] - 1, commas[4 - 1] + 1);
    mu = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[5] - commas[5 - 1] - 1, commas[5 - 1] + 1);
    a = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, line.size(), commas[6 - 1] + 1);
    S = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.clear();

    getline(ifs, line);
    tmp = int(line.size());
    for (int j = 0; j < tmp; j++) {
        if (line.compare(j, 1, ",") == 0) {
            commas.push_back(j);
            contc++;
        }
    }
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[0], 0);
    center.push_back(atof(tmpme2));
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[1] - commas[1 - 1] - 1, commas[1 - 1] + 1);
    center.push_back(atof(tmpme2));
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[2] - commas[2 - 1] - 1, commas[2 - 1] + 1);
    center.push_back(atof(tmpme2));
    memset(tmpme2, 0, maxsize);
    line.clear();
#ifdef PARALLEL
    vector<int>::iterator it;
    int rank;
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    while (getline(ifs, line)) {
        if (line.compare(0, 1, "*") == 0) {
            ifs.seekg(oldpos);
            break;
        } else {
            //PARALLEL
            if (epart.size() != 0) {
#ifdef PARALLEL
                int size = int(line.size());
                line.copy(tmpme2, size, 0);
                vector<classElements *> subElements;
                int mainEle = atoi(tmpme2) - 1;
                memset(tmpme2, 0, maxsize);
                if (epart[mainEle] == rank) {
                    it = find(elemLocal.begin(), elemLocal.end(), mainEle);
                    mainEle = it - elemLocal.begin();
                    elements[mainEle]->getSubElement(elemsPress.size(), S, subElements, nodes, mainEle,
                                                     surface_order); //is me a unique identifier? dont think it is being used or matters
                                                     for (int i = 0; i < subElements.size(); i++) {
                                                         elemsPress.push_back(subElements[i]);
                                                     }

                }
#endif
                //SERIAL
            } else {
                int size = int(line.size());
                line.copy(tmpme2, size, 0);
                vector<classElements *> subElements;
                int mainEle = atoi(tmpme2) - 1;
                memset(tmpme2, 0, maxsize);
                elements[mainEle]->getSubElement(elemsPress.size(), S, subElements, nodes, mainEle,
                                                 surface_order); //is me a unique identifier? dont think it is being used or matters
                                                 for (int i = 0; i < subElements.size(); i++) {
                                                     elemsPress.push_back(subElements[i]);
                                                 }
            }
        }
        oldpos = ifs.tellg();
    }
    nodForces = new classHertzianRamp(type, elemsPress, P, T, mu, a, center, t0, tf);
    NeumannBCs.push_back(nodForces);

}
/*! \brief Load the extra Neumann BCs from the input file
  @param[in] elements Array with all elements (or background cells) in the domain
  @param[in] nodes Array of nodes
  @param[inout] NeumannBCs Array with the Neumann boundary conditions
  @param[in] ndim Dimension of the problem
  @param[in] ifs Pointer to the input file (.inp)
  @param[in] npart Processor ownership mapping of the nodes in the whole model npart[node i]=rank (empty for sequential simulations)
  @param[in] nodLocal Mapping of local nodes to global nodes numbering nodlocal[local node index]= global node index (empty for sequential simulations)
*/
extern void currentInstBCManagement(vector<classElements *> &elements, vector<classNodes *> &nodes,
                                    vector<classNeumannBCs *> &NeumannBCs, int ndim, ifstream &ifs, vector<int> &epart,
                                    vector<int> &elemLocal) {

    classNeumannBCs *nodBCs;
    vector<classElements *> elemsCurrent;
    string type = "CURRENT INST";

    vector<int> commas;
    int maxsize = 100;
    char tmpme2[maxsize], tempU[maxsize];
    memset(tmpme2, 0, maxsize);
    memset(tempU, 0, maxsize);
    double tmp = 0;
    string line;
    streampos oldpos;

    getline(ifs, line);
    tmp = int(line.size());
    int contc = 0;
    for (int j = 0; j < tmp; j++) {
        if (line.compare(j, 1, ",") == 0) {
            commas.push_back(j);
            contc++;
        }
    }

    // Parameters
    double t0 = 0;
    double tf = 0;
    double INeumann = 0;

    // Reading the parameters
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[0], 0);
    t0 = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[1] - commas[1 - 1] - 1, commas[1 - 1] + 1);
    tf = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, line.size(), commas[2 - 1] + 1);
    INeumann = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.clear();

#ifdef PARALLEL
    vector<int>::iterator it;
    int rank;
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    // Getting the elements
    while (getline(ifs, line)) {
        if (line.compare(0, 1, "*") == 0) {
            ifs.seekg(oldpos);
            break;
        } else {
            int size = int(line.size());
            line.copy(tmpme2, size, 0);
            int Element = atoi(tmpme2) - 1;
            memset(tmpme2, 0, maxsize);
#ifdef PARALLEL
            if (epart[Element] == rank) {
                it = find(elemLocal.begin(), elemLocal.end(), Element);
                Element = it - elemLocal.begin();
                elemsCurrent.push_back(elements[Element]);
            }
#else
            elemsCurrent.push_back(elements[Element]);
#endif
        }
        oldpos = ifs.tellg();
    }

    // Creating the classCurrentInst object
    nodBCs = new classCurrentInst(type, elemsCurrent, INeumann, t0, tf);
    NeumannBCs.push_back(nodBCs);

}

/*! \brief Load the extra Neumann BCs (Volumetric heat flux) from the input file
  @param[in] elements Array with all elements (or background cells) in the domain
  @param[in] nodes Array of nodes
  @param[inout] NeumannBCs Array with the Neumann boundary conditions
  @param[in] ndim Dimension of the problem
  @param[in] ifs Pointer to the input file (.inp)
  @param[in] npart Processor ownership mapping of the nodes in the whole model npart[node i]=rank (empty for sequential simulations)
  @param[in] nodLocal Mapping of local nodes to global nodes numbering nodlocal[local node index]= global node index (empty for sequential simulations)
*/
extern void volHeatBCManagement(vector<classElements *> &elements, vector<classNodes *> &nodes,
                                vector<classNeumannBCs *> &NeumannBCs, int ndim, ifstream &ifs, vector<int> &epart,
                                vector<int> &elemLocal) {

    classNeumannBCs *nodBCs;
    vector<classElements *> elemsVolHeat;
    string type = "VOLUMETRIC HEAT FLUX";

    vector<int> commas;
    int maxsize = 100;
    char tmpme2[maxsize], tempU[maxsize];
    memset(tmpme2, 0, maxsize);
    memset(tempU, 0, maxsize);
    double tmp = 0;
    string line;
    streampos oldpos;

    getline(ifs, line);
    tmp = int(line.size());
    int contc = 0;
    for (int j = 0; j < tmp; j++) {
        if (line.compare(j, 1, ",") == 0) {
            commas.push_back(j);
            contc++;
        }
    }

    // Parameters
    double t0 = 0;
    double tf = 0;
    double r = 0;

    // Reading the parameters
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[0], 0);
    t0 = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[1] - commas[1 - 1] - 1, commas[1 - 1] + 1);
    tf = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, line.size(), commas[2 - 1] + 1);
    r = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.clear();

#ifdef PARALLEL
    vector<int>::iterator it;
    int rank;
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    // Getting the elements
    while (getline(ifs, line)) {
        if (line.compare(0, 1, "*") == 0) {
            ifs.seekg(oldpos);
            break;
        } else {
            int size = int(line.size());
            line.copy(tmpme2, size, 0);
            int Element = atoi(tmpme2) - 1;
            memset(tmpme2, 0, maxsize);
#ifdef PARALLEL
            if (epart[Element] == rank) {
                it = find(elemLocal.begin(), elemLocal.end(), Element);
                Element = it - elemLocal.begin();
                elemsVolHeat.push_back(elements[Element]);
            }
#else
            elemsVolHeat.push_back(elements[Element]);
#endif
        }
        oldpos = ifs.tellg();
    }

    // Creating the classVolHeat object
    nodBCs = new classVolHeatFlux(type, elemsVolHeat, r, t0, tf);
    NeumannBCs.push_back(nodBCs);
    //delete nodBCs;
}

/*! \brief Load the extra Neumann BCs (Volumetric heat flux) from the input file
  @param[in] elements Array with all elements (or background cells) in the domain
  @param[in] nodes Array of nodes
  @param[inout] NeumannBCs Array with the Neumann boundary conditions
  @param[in] ndim Dimension of the problem
  @param[in] ifs Pointer to the input file (.inp)
  @param[in] npart Processor ownership mapping of the nodes in the whole model npart[node i]=rank (empty for sequential simulations)
  @param[in] nodLocal Mapping of local nodes to global nodes numbering nodlocal[local node index]= global node index (empty for sequential simulations)
*/
extern void gaussianVolHeatBCManagement(vector<classElements *> &elements, vector<classNodes *> &nodes,
                                vector<classNeumannBCs *> &NeumannBCs, int ndim, ifstream &ifs, vector<int> &epart,
                                vector<int> &elemLocal) {

    string type = "GAUSSIAN VOLUMETRIC HEAT FLUX";
    
    vector<int> commas;
    int maxsize = 100;
    char tmpme2[maxsize], tempU[maxsize];
    memset(tmpme2, 0, maxsize);
    memset(tempU, 0, maxsize);
    double tmp = 0;
    string line;
    streampos oldpos;

    getline(ifs, line);
    tmp = int(line.size());
    int contc = 0;
    for (int j = 0; j < tmp; j++) {
        if (line.compare(j, 1, ",") == 0) {
            commas.push_back(j);
            contc++;
        }
    }

    // Parameters
    double t0 = 0;
    double tf = 0;
    double qmax = 0;
    double x0(0.), y0(0.), z0(0.);
    double vx(0), vy(0), vz(0);
    double a(0.), b(0.), c(0.);

    // Reading the parameters
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[0], 0);
    t0 = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[1] - commas[1 - 1] - 1, commas[1 - 1] + 1);
    tf = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[2] - commas[2 - 1] - 1, commas[2 - 1] + 1);
    qmax = atof(tmpme2);
    //
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[3] - commas[3 - 1] - 1, commas[3 - 1] + 1);
    x0 = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[4] - commas[4 - 1] - 1, commas[4 - 1] + 1);
    y0 = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[5] - commas[5 - 1] - 1, commas[5 - 1] + 1);
    z0 = atof(tmpme2);
    //
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[6] - commas[6 - 1] - 1, commas[6 - 1] + 1);
    vx = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[7] - commas[7 - 1] - 1, commas[7 - 1] + 1);
    vy = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[8] - commas[8 - 1] - 1, commas[8 - 1] + 1);
    vz = atof(tmpme2);
    //
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[9] - commas[9 - 1] - 1, commas[9 - 1] + 1);
    a = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[10] - commas[10 - 1] - 1, commas[10 - 1] + 1);
    b = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, line.size(), commas[11 - 1] + 1);
    c = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.clear();
    // Creating the classVolHeat object
    classNeumannBCs* BC = new classVolHeatFluxGaussian(type, qmax, x0, y0, z0, vx, vy, vz, a, b, c , t0, tf);
    NeumannBCs.push_back(BC);
}


extern void stochasticHertzianRampBCManagement(vector<classElements *> &elements, vector<classNodes *> &nodes,
                                     vector<classNeumannBCs *> &NeumannBCs, int ndim, ifstream &ifs, vector<int> &epart,
                                     vector<int> &elemLocal, int surface_order) {

    classNeumannBCs *nodForces;
    string type = "STOCHASTIC HERTZIAN RAMP";
    bool flagHaar = false;
    string distribution= "";
    int resolution = 0;
    vector<double> tempxyz;
    vector<int> commas;
    int maxsize = 100;
    char tmpme2[maxsize], tempU[maxsize];
    memset(tmpme2, 0, maxsize);
    memset(tempU, 0, maxsize);
    double tmp = 0;
    string line;
    double P = 0;
    double T = 0;
    double mu;
    double mu_var;
    double a;
    int S = 0;
    vector<double> center;
    vector<classElements *> elemsPress;
    streampos oldpos;

    getline(ifs, line);
    tmp = int(line.size());
    int contc = 0;
    for (int j = 0; j < tmp; j++) {
        if (line.compare(j, 1, ",") == 0) {
            commas.push_back(j);
            contc++;
        }
    }

    double t0 = 0;
    double tf = 0;
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[0], 0);
    t0 = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[1] - commas[1 - 1] - 1, commas[1 - 1] + 1);
    tf = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[2] - commas[2 - 1] - 1, commas[2 - 1] + 1);
    P = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[3] - commas[3 - 1] - 1, commas[3 - 1] + 1);
    T = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[4] - commas[4 - 1] - 1, commas[4 - 1] + 1);
    mu = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[5] - commas[5 - 1] - 1, commas[5 - 1] + 1);
    mu_var = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[6] - commas[6 - 1] - 1, commas[6 - 1] + 1);
    a = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, line.size(), commas[7 - 1] + 1);
    S = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.clear();

    getline(ifs, line);
    tmp = int(line.size());
    for (int j = 0; j < tmp; j++) {
        if (line.compare(j, 1, ",") == 0) {
            commas.push_back(j);
            contc++;
        }
    }
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[0], 0);
    center.push_back(atof(tmpme2));
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[1] - commas[1 - 1] - 1, commas[1 - 1] + 1);
    center.push_back(atof(tmpme2));
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[2] - commas[2 - 1] - 1, commas[2 - 1] + 1);
    center.push_back(atof(tmpme2));
    memset(tmpme2, 0, maxsize);
    line.clear();
#ifdef PARALLEL
    vector<int>::iterator it;
    int rank;
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    while (getline(ifs, line)) {
        if (line.compare(0, 1, "*") == 0) {
            ifs.seekg(oldpos);
            break;
        } else {
            //PARALLEL
            if (epart.size() != 0) {
#ifdef PARALLEL
                int size = int(line.size());
                line.copy(tmpme2, size, 0);
                vector<classElements *> subElements;
                int mainEle = atoi(tmpme2) - 1;
                memset(tmpme2, 0, maxsize);
                if (epart[mainEle] == rank) {
                    it = find(elemLocal.begin(), elemLocal.end(), mainEle);
                    mainEle = it - elemLocal.begin();
                    elements[mainEle]->getSubElement(elemsPress.size(), S, subElements, nodes, mainEle,
                                                     surface_order); //is me a unique identifier? dont think it is being used or matters
                                                     for (int i = 0; i < subElements.size(); i++) {
                                                         elemsPress.push_back(subElements[i]);
                                                     }

                }
#endif
                //SERIAL
            } else {
                int size = int(line.size());
                line.copy(tmpme2, size, 0);
                vector<classElements *> subElements;
                int mainEle = atoi(tmpme2) - 1;
                memset(tmpme2, 0, maxsize);
                elements[mainEle]->getSubElement(elemsPress.size(), S, subElements, nodes, mainEle,
                                                 surface_order); //is me a unique identifier? dont think it is being used or matters
                                                 for (int i = 0; i < subElements.size(); i++) {
                                                     elemsPress.push_back(subElements[i]);
                                                 }
            }
        }
        oldpos = ifs.tellg();
    }
    vector<constitutiveModels *> elementsCons = elements[0]->getConsModel();
    int nbr_dof;
    nbr_dof = (elementsCons[0]->getNbrDofConsMod());
    int nstoch;
    nstoch = nbr_dof/ndim;
    vector<double> C(nstoch*nstoch*nstoch);
    const classStochasticHyperElasticStVenantKirchhoff* firstLaw = dynamic_cast<const classStochasticHyperElasticStVenantKirchhoff*>(elementsCons[0]);
    if (firstLaw == NULL)
    {
        ERROR("classStochasticHyperElasticStVenantKirchhoff must be used");
        exit(-1);
    }
    C = firstLaw->getC_Stoch();
    while (getline(ifs, line)) {
        if(line.compare(0,6,"*ORDER") == 0) {
            getline(ifs, line);
            int size = int(line.size());
            line.copy(tmpme2, size, 0);
            resolution = atof(tmpme2);
            memset(tmpme2, 0, maxsize);
            line.clear();
        }
        if (line.compare(0, 4, "Haar") == 0) {
            flagHaar = true;
        }
        if (line.compare(0, 13, "*DISTRIBUTION") == 0){
            getline(ifs, line);        
                    for (int i = 0; i < line.size(); i++) {
                        if (line[i] != ' ')//remove all whitespaces
                        {
                            distribution = distribution + line[i];
                        }
                    }
        }
    }
    nodForces = new classStochasticHertzianRamp(type, elemsPress, P, T, mu, mu_var, a, center, t0, tf, nstoch, C, flagHaar, resolution, distribution);
    NeumannBCs.push_back(nodForces);

}



/*! \brief Load the extra Neumann BCs (Surface heat flux) from the input file
  @param[in] elements Array with all elements (or background cells) in the domain
  @param[in] nodes Array of nodes
  @param[inout] NeumannBCs Array with the Neumann boundary conditions
  @param[in] ndim Dimension of the problem
  @param[in] ifs Pointer to the input file (.inp)
  @param[in] npart Processor ownership mapping of the nodes in the whole model npart[node i]=rank (empty for sequential simulations)
  @param[in] nodLocal Mapping of local nodes to global nodes numbering nodlocal[local node index]= global node index (empty for sequential simulations)
*/
extern void heatInstBCManagement(vector<classElements *> &elements, vector<classNodes *> &nodes,
                                 vector<classNeumannBCs *> &NeumannBCs, int ndim, ifstream &ifs, vector<int> &epart,
                                 vector<int> &elemLocal, int surface_order) {

    classNeumannBCs *nodForces;
    string type = "HEAT INST";
    vector<double> tempxyz;
    vector<int> commas;
    int maxsize = 100;
    char tmpme2[maxsize], tempU[maxsize];
    memset(tmpme2, 0, maxsize);
    memset(tempU, 0, maxsize);
    double tmp = 0;
    string line;
    vector<double> heatflux;
    int S = 0;
    vector<classElements *> elemsHeat;
    streampos oldpos;

    getline(ifs, line);
    tmp = int(line.size());
    int contc = 0;
    for (int j = 0; j < tmp; j++) {
        if (line.compare(j, 1, ",") == 0) {
            commas.push_back(j);
            contc++;
        }
    }

    double t0 = 0;
    double tf = 0;
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[0], 0);
    t0 = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[1] - commas[1 - 1] - 1, commas[1 - 1] + 1);
    tf = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[2] - commas[2 - 1] - 1, commas[2 - 1] + 1);
    heatflux.push_back(atof(tmpme2));
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[3] - commas[3 - 1] - 1, commas[3 - 1] + 1);
    heatflux.push_back(atof(tmpme2));
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[4] - commas[4 - 1] - 1, commas[4 - 1] + 1);
    heatflux.push_back(atof(tmpme2));
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, line.size(), commas[5 - 1] + 1);
    S = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.clear();
#ifdef PARALLEL
    vector<int>::iterator it;
    int rank;
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    while (getline(ifs, line)) {
        if (line.compare(0, 1, "*") == 0) {
            ifs.seekg(oldpos);
            break;
        } else {
            //PARALLEL
            if (epart.size() != 0) {
#ifdef PARALLEL
                int size = int(line.size());
                line.copy(tmpme2, size, 0);
                vector<classElements *> subElements;
                int mainEle = atoi(tmpme2) - 1;
                memset(tmpme2, 0, maxsize);
                if (epart[mainEle] == rank) {
                    it = find(elemLocal.begin(), elemLocal.end(), mainEle);
                    mainEle = it - elemLocal.begin();
                    elements[mainEle]->getSubElement(elemsHeat.size(), S, subElements, nodes, mainEle,
                                                     surface_order); //is me a unique identifier? dont think it is being used or matters
                    for (int i = 0; i < subElements.size(); i++) {
                        elemsHeat.push_back(subElements[i]);
                    }
                }
#endif
                //SERIAL
            } else {
                int size = int(line.size());
                line.copy(tmpme2, size, 0);
                vector<classElements *> subElements;
                int mainEle = atoi(tmpme2) - 1;
                memset(tmpme2, 0, maxsize);
                elements[mainEle]->getSubElement(elemsHeat.size(), S, subElements, nodes, mainEle,
                                                 surface_order); //is me a unique identifier? dont think it is being used or matters
                for (int i = 0; i < subElements.size(); i++) {
                    elemsHeat.push_back(subElements[i]);
                }
            }
        }
        oldpos = ifs.tellg();
    }
    nodForces = new classSurfaceHeatFluxInst(type, elemsHeat, heatflux, t0, tf);
    NeumannBCs.push_back(nodForces);

}

/*! \brief Load the extra Convection BCs from the input file
  @param[in] elements Array with all elements (or background cells) in the domain
  @param[in] nodes Array of nodes
  @param[inout] NeumannBCs Array with the Neumann boundary conditions
  @param[in] ndim Dimension of the problem
  @param[in] ifs Pointer to the input file (.inp)
  @param[in] npart Processor ownership mapping of the nodes in the whole model npart[node i]=rank (empty for sequential simulations)
  @param[in] nodLocal Mapping of local nodes to global nodes numbering nodlocal[local node index]= global node index (empty for sequential simulations)
*/
extern void convectionBCManagement(vector<classElements *> &elements, vector<classNodes *> &nodes,
                                   vector<classNeumannBCs *> &NeumannBCs, int ndim, ifstream &ifs, vector<int> &epart,
                                   vector<int> &elemLocal, int surface_order) {

    classNeumannBCs *nodForces;
    string type = "CONVECTION";
    vector<double> tempxyz;
    vector<int> commas;
    int maxsize = 100;
    char tmpme2[maxsize], tempU[maxsize];
    memset(tmpme2, 0, maxsize);
    memset(tempU, 0, maxsize);
    double tmp = 0;
    string line;
    int S = 0;
    vector<classElements *> elemsConvection;
    streampos oldpos;

    getline(ifs, line);
    tmp = int(line.size());
    int contc = 0;
    for (int j = 0; j < tmp; j++) {
        if (line.compare(j, 1, ",") == 0) {
            commas.push_back(j);
            contc++;
        }
    }

    double t0 = 0;
    double tf = 0;
    double filmCoefficient;
    double sinkTemp;
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[0], 0);
    t0 = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[1] - commas[1 - 1] - 1, commas[1 - 1] + 1);
    tf = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[2] - commas[2 - 1] - 1, commas[2 - 1] + 1);
    filmCoefficient = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[3] - commas[3 - 1] - 1, commas[3 - 1] + 1);
    sinkTemp = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, line.size(), commas[4 - 1] + 1);
    S = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.clear();
#ifdef PARALLEL
    vector<int>::iterator it;
    int rank;
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    while (getline(ifs, line)) {
        if (line.compare(0, 1, "*") == 0) {
            ifs.seekg(oldpos);
            break;
        } else {
            //PARALLEL
            if (epart.size() != 0) {
#ifdef PARALLEL
                int size = int(line.size());
                line.copy(tmpme2, size, 0);
                vector<classElements *> subElements;
                int mainEle = atoi(tmpme2) - 1;
                memset(tmpme2, 0, maxsize);
                if (epart[mainEle] == rank) {
                    it = find(elemLocal.begin(), elemLocal.end(), mainEle);
                    mainEle = it - elemLocal.begin();
                    elements[mainEle]->getSubElement(elemsConvection.size(), S, subElements, nodes, mainEle,
                                                     surface_order); //is me a unique identifier? dont think it is being used or matters
                    for (int i = 0; i < subElements.size(); i++) {
                        elemsConvection.push_back(subElements[i]);
                    }
                }
#endif
                //SERIAL
            } else {
                int size = int(line.size());
                line.copy(tmpme2, size, 0);
                vector<classElements *> subElements;
                int mainEle = atoi(tmpme2) - 1;
                memset(tmpme2, 0, maxsize);
                elements[mainEle]->getSubElement(elemsConvection.size(), S, subElements, nodes, mainEle,
                                                 surface_order); //is me a unique identifier? dont think it is being used or matters
                for (int i = 0; i < subElements.size(); i++) {
                    elemsConvection.push_back(subElements[i]);
                }
            }
        }
        oldpos = ifs.tellg();
    }
    nodForces = new classConvection(type, elemsConvection, filmCoefficient, sinkTemp, t0, tf);
    NeumannBCs.push_back(nodForces);

}

/*! \brief Load the extra Radiation BCs from the input file
  @param[in] elements Array with all elements (or background cells) in the domain
  @param[in] nodes Array of nodes
  @param[inout] NeumannBCs Array with the Neumann boundary conditions
  @param[in] ndim Dimension of the problem
  @param[in] ifs Pointer to the input file (.inp)
  @param[in] npart Processor ownership mapping of the nodes in the whole model npart[node i]=rank (empty for sequential simulations)
  @param[in] nodLocal Mapping of local nodes to global nodes numbering nodlocal[local node index]= global node index (empty for sequential simulations)
*/
extern void radiationBCManagement(vector<classElements *> &elements, vector<classNodes *> &nodes,
                                  vector<classNeumannBCs *> &NeumannBCs, int ndim, ifstream &ifs, vector<int> &epart,
                                  vector<int> &elemLocal, int surface_order) {

    classNeumannBCs *nodForces;
    string type = "RADIATION";
    vector<double> tempxyz;
    vector<int> commas;
    int maxsize = 100;
    char tmpme2[maxsize], tempU[maxsize];
    memset(tmpme2, 0, maxsize);
    memset(tempU, 0, maxsize);
    double tmp = 0;
    string line;
    int S = 0;
    vector<classElements *> elemsConvection;
    streampos oldpos;

    getline(ifs, line);
    tmp = int(line.size());
    int contc = 0;
    for (int j = 0; j < tmp; j++) {
        if (line.compare(j, 1, ",") == 0) {
            commas.push_back(j);
            contc++;
        }
    }

    double t0 = 0;
    double tf = 0;
    double radiationConstant;
    double sinkTemp;
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[0], 0);
    t0 = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[1] - commas[1 - 1] - 1, commas[1 - 1] + 1);
    tf = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[2] - commas[2 - 1] - 1, commas[2 - 1] + 1);
    radiationConstant = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, commas[3] - commas[3 - 1] - 1, commas[3 - 1] + 1);
    sinkTemp = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.copy(tmpme2, line.size(), commas[4 - 1] + 1);
    S = atof(tmpme2);
    memset(tmpme2, 0, maxsize);
    line.clear();
#ifdef PARALLEL
    vector<int>::iterator it;
    int rank;
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    while (getline(ifs, line)) {
        if (line.compare(0, 1, "*") == 0) {
            ifs.seekg(oldpos);
            break;
        } else {
            //PARALLEL
            if (epart.size() != 0) {
#ifdef PARALLEL
                int size = int(line.size());
                line.copy(tmpme2, size, 0);
                vector<classElements *> subElements;
                int mainEle = atoi(tmpme2) - 1;
                memset(tmpme2, 0, maxsize);
                if (epart[mainEle] == rank) {
                    it = find(elemLocal.begin(), elemLocal.end(), mainEle);
                    mainEle = it - elemLocal.begin();
                    elements[mainEle]->getSubElement(elemsConvection.size(), S, subElements, nodes, mainEle,
                                                     surface_order); //is me a unique identifier? dont think it is being used or matters
                    for (int i = 0; i < subElements.size(); i++) {
                        elemsConvection.push_back(subElements[i]);
                    }
                }
#endif
                //SERIAL
            } else {
                int size = int(line.size());
                line.copy(tmpme2, size, 0);
                vector<classElements *> subElements;
                int mainEle = atoi(tmpme2) - 1;
                memset(tmpme2, 0, maxsize);
                elements[mainEle]->getSubElement(elemsConvection.size(), S, subElements, nodes, mainEle,
                                                 surface_order); //is me a unique identifier? dont think it is being used or matters
                for (int i = 0; i < subElements.size(); i++) {
                    elemsConvection.push_back(subElements[i]);
                }
            }
        }
        oldpos = ifs.tellg();
    }
    nodForces = new classRadiation(type, elemsConvection, radiationConstant, sinkTemp, t0, tf);
    NeumannBCs.push_back(nodForces);

}

/*! \brief It reads all Neumann BCs
  @param[in] elements Array with all elements (or background cells) in the domain
  @param[in] nodes Array with all nodes in the domain
  @param[out] NeumannBCs Array with the Neumann boundary conditions
  @param[in] ndim Dimension of the problem
  @param[in] ifs Pointer to the input file (.inp)
  @param[in] npart Processor ownership mapping of the nodes in the whole model npart[node i]=rank (empty for sequential simulations)
  @param[in] nodLocal Mapping of local nodes to global nodes numbering nodlocal[local node index]= global node index (empty for sequential simulations)
*/
extern void NeumannCommonManagement(vector<classElements *> &elements, vector<classNodes *> &nodes,
                                    vector<classNeumannBCs *> &NeumannBCs, int ndim, ifstream &ifs, vector<int> &npart,
                                    vector<int> &nodLocal, vector<int> &epart, vector<int> &elemLocal,
                                    int surface_order) {
    string line;
    while (getline(ifs, line)) {
        if (line.compare(0, 9, "*BOUNDARY") == 0) {
            int tmp2 = int(line.size());
            if (line.compare(16, tmp2, "FORCE") == 0) {
                forceBCManagement(nodes, NeumannBCs, ndim, ifs, npart, nodLocal);
            } else if (line.compare(16, 29, "TRACTION RAMP") == 0) {
                tractionRampBCManagement(elements, nodes, NeumannBCs, ndim, ifs, epart, elemLocal, surface_order);
            } else if (line.compare(16, 29, "TRACTION INST") == 0) {
                tractionInstBCManagement(elements, nodes, NeumannBCs, ndim, ifs, epart, elemLocal, surface_order);
            } else if (line.compare(16, 29, "PRESSURE RAMP") == 0) {
                pressureRampBCManagement(elements, nodes, NeumannBCs, ndim, ifs, epart, elemLocal, surface_order);
            } else if (line.compare(16, 29, "PRESSURE INST") == 0) {
                pressureInstBCManagement(elements, nodes, NeumannBCs, ndim, ifs, epart, elemLocal, surface_order);
            }
        }
        if (line.compare(0, 8, "*GRAVITY") == 0) {
            int maxsize = 100;
            double g;
            char tempNode[maxsize];
            memset(tempNode, 0, maxsize);
            classNeumannBCs *nodBCs;
            string type = "GRAVITY";
            getline(ifs, line);
            line.copy(tempNode, line.size(), 0);
            g = atof(tempNode);
            nodBCs = new classGravity(type, g);
            NeumannBCs.push_back(nodBCs);

        }
        if (line.compare(0, 17, "*HERTZIAN CONTACT") == 0) {
            hertzianRampBCManagement(elements, nodes, NeumannBCs, ndim, ifs, epart, elemLocal, surface_order);
        }
        if (line.compare(0, 28, "*STOCHASTIC HERTZIAN CONTACT") == 0) {
            stochasticHertzianRampBCManagement(elements, nodes, NeumannBCs, ndim, ifs, epart, elemLocal, surface_order);
        }
        if (line.compare(0, 14, "*FLUX EXTRADOF") == 0) {
            if (line.compare(21, 33, "CURRENT INST") == 0) {
                currentInstBCManagement(elements, nodes, NeumannBCs, ndim, ifs, epart, elemLocal);
            } else if (line.compare(21, 30, "HEAT INST") == 0) {
                heatInstBCManagement(elements, nodes, NeumannBCs, ndim, ifs, epart, elemLocal, surface_order);
            } else if (line.compare(21, 31, "CONVECTION") == 0) {
                convectionBCManagement(elements, nodes, NeumannBCs, ndim, ifs, epart, elemLocal, surface_order);
            } else if (line.compare(21, 30, "RADIATION") == 0) {
                radiationBCManagement(elements, nodes, NeumannBCs, ndim, ifs, epart, elemLocal, surface_order);
            } else if (line.compare(21, 41, "VOLUMETRIC HEAT FLUX") == 0) {
                volHeatBCManagement(elements, nodes, NeumannBCs, ndim, ifs, epart, elemLocal);
            } else if (line.compare(21, 45, "GAUSSIAN VOLUMETRIC HEAT FLUX") == 0) {
                gaussianVolHeatBCManagement(elements, nodes, NeumannBCs, ndim, ifs, epart, elemLocal);
            }
        }
        if (line.compare(0, 4, "*FSI") == 0) {
            FSIManagement(nodes,NeumannBCs, ifs, npart,nodLocal);
        }
    }
}


/*! \brief It reads the the number of outputs that the user wants
  @param[out] outputs Number of outputs
  @param[in] ifs Pointer to the input file (.inp)
*/
extern void outManagement(classOutputs *&outputs, ifstream &ifs) {

    streampos oldpos;
    oldpos = ifs.tellg();

    string line, secondLine, thirdLine;
    int maxsize = 100;
    char tmpme2[maxsize];
    vector<int> commas;
    vector<int> outVars;

    getline(ifs, secondLine);
    outputs->setNumFiles(atoi(secondLine.c_str()));
    getline(ifs, thirdLine);
    if (thirdLine.compare(0, 1, "*") != 0) {
        int tmp = int(thirdLine.size());
        if (thirdLine.compare(0, 3, "All") == 0) {
            memset(tmpme2, 0, maxsize);
            thirdLine.copy(tmpme2, tmp, 4);
            int max = atoi(tmpme2);
            for (int j = 0; j < max; j++) {
                outVars.push_back(j);
            }
            outputs->setOutVars(outVars);
        } else {
            int contc = 0;
            commas.clear();
            for (int j = 0; j < tmp; j++) {
                if (thirdLine.compare(j, 1, ",") == 0) {
                    commas.push_back(j);
                    contc++;
                }
            }
            if (commas.size() == 0) {
                outVars.push_back(atoi(thirdLine.c_str()));
                outputs->setOutVars(outVars);
            } else {
                memset(tmpme2, 0, maxsize);
                thirdLine.copy(tmpme2, commas[0], 0);
                int max = atoi(tmpme2); // Maximum number of internal variables
                for (int k = 1; k < commas.size(); k++) {
                    memset(tmpme2, 0, maxsize);
                    thirdLine.copy(tmpme2, commas[k] - commas[k - 1] - 1, commas[k - 1] + 1);
                    outVars.push_back(atoi(tmpme2) - 1);
                }
                thirdLine.copy(tmpme2, tmp - commas[commas.size() - 1] - 1, commas[commas.size() - 1] + 1);
                outVars.push_back(atoi(tmpme2) - 1);
                outputs->setOutVars(outVars);
                for (int i = 0; i < outVars.size(); i++) {
                    if (outVars[i] > max) {
                        ERROR("The maximum number of internal variables in the whole domain is %i, and you asked for the position %i to be in the outputs files",
                              max, outVars[i]);
                        exit(0);
                    }
                }
            }
        }
    }
    ifs.seekg(oldpos);
}

extern void forceOutManagement(classPrintForces *&forceOutputs, ifstream &ifs) {
    int maxsize = 100;
    char tmpme2[maxsize], tempDirection[maxsize];
    memset(tmpme2, 0, maxsize);
    memset(tempDirection, 0, maxsize);
    string line, templine2;
    int tmp = 0;
    char tempNode[maxsize];
    memset(tempNode, 0, maxsize);
    vector<int> commas;
    int tempme = 0;
    streampos oldpos;

    vector<int> nodesToWrite;
    vector<int> directionsToWrite;

    templine2.clear();
    getline(ifs, templine2);
    int size = int(templine2.size());
    templine2.copy(tempNode, size, 0);
    double frequency = atof(tempNode);
    memset(tempNode, 0, maxsize);
    templine2.clear();
    forceOutputs->setPeriod(1 / frequency);


    while (getline(ifs, line)) {
        if (line.compare(0, 1, "*") == 0) {
            ifs.seekg(oldpos);
            break;
        } else {
            tmp = int(line.size());
            int contc = 0;
            commas.clear();
            for (int j = 0; j < line.size(); j++) {
                if (line.compare(j, 1, ",") == 0) {
                    commas.push_back(j);
                    contc++;
                }
            }
            line.copy(tmpme2, commas[0], 0);
            // What direction
            // 1-x, 2-y, 3-z
            line.copy(tempDirection, commas[1] - commas[1 - 1] - 1, commas[1 - 1] + 1);
            tempme = atoi(tmpme2) - 1;
            int direction = atoi(tempDirection);

            nodesToWrite.push_back(tempme);
            directionsToWrite.push_back(direction);
            memset(tmpme2, 0, maxsize);
            memset(tempDirection, 0, maxsize);
            commas.clear();
        }
        oldpos = ifs.tellg();
    }
    forceOutputs->setDirectionsToWrite(directionsToWrite);
    forceOutputs->setNodesToWrite(nodesToWrite);

}

/*! \brief It reads the FSI nodes in the case an FSI solution is desired
  @param[inout] NeumannBCs Array with the Neumann boundary conditions
  @param[in] ifs Pointer to the input file (.inp)
  @param[in] npart Processor ownership of the nodes in the whole model npartShare[node i]=1 if the nodes is mine or
  shared with other processor, 0 otherwise (empty for sequential simulations)
*/
extern void FSIManagement(vector<classNodes *> &nodes, vector<classNeumannBCs *> &NeumannBCs, ifstream &ifs, vector<int> npart, vector<int> &nodLocal) {

    int maxsize = 100;
    string type = "FSI";
    classNeumannBCs *FSIBC;
    vector<int> boundary_nodes;
    vector<double> tempxyz;
    char tempNode[maxsize];
    string line;
    int tmp = 0, size = 0;
    vector<vector<double> > boundary_coordinates;

#ifdef PARALLEL
    vector<int>::iterator it;
    int rank;
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    while (getline(ifs, line)) {
        if (line.compare(0, 1, "*") == 0) {
            break;
        } else {
            size = int(line.size());
            line.copy(tempNode, size, 0);
            tmp = atoi(tempNode) - 1;
            
#ifdef PARALLEL
            if (npart[tmp] == rank) {
                it = find(nodLocal.begin(), nodLocal.end(), tmp);
                if (it == nodLocal.end())
                {
                  ERROR("node %d is not found",tmp);
                }
                else
                {
                  tempxyz = nodes[(it - nodLocal.begin())]->getXYZ();
                  boundary_coordinates.push_back(tempxyz);
                  boundary_nodes.push_back(it - nodLocal.begin());
                }
            }
#else
            tempxyz = nodes[tmp]->getXYZ();
            boundary_coordinates.push_back(tempxyz);
            boundary_nodes.push_back(tmp);        
#endif
            memset(tempNode, 0, maxsize);
        }
    }
    FSIBC = new classFSI(type, boundary_nodes, boundary_coordinates);
    NeumannBCs.push_back(FSIBC);
}


extern void NodeDataExtractionManagement(vector<classNodes *> &nodes, ifstream &ifs, vector<int> npart, vector<int> &nodLocal)
{
    int maxsize = 100;
    char tempNode[maxsize];
    int tmp = 0, size = 0;
    vector<int> commas;
    string line;
    getline(ifs, line);
    tmp = int(line.size());
    int contc = 0;
    for (int j = 0; j < tmp; j++) {
        if (line.compare(j, 1, ",") == 0) {
            commas.push_back(j);
            contc++;
        }
    }
    vector<int> allNodes;
    string name, dataType, operation;
    int comp;
    
    memset(tempNode, 0, maxsize);
    line.copy(tempNode, commas[0], 0);
    name = tempNode;
    memset(tempNode, 0, maxsize);
    line.copy(tempNode, commas[1] - commas[1 - 1] - 1, commas[1 - 1] + 1);
    dataType = tempNode;
    memset(tempNode, 0, maxsize);
    line.copy(tempNode, commas[2] - commas[2 - 1] - 1, commas[2 - 1] + 1);
    operation = tempNode;
    memset(tempNode, 0, maxsize);
    line.copy(tempNode, line.size(), commas[3 - 1] + 1);
    comp = atof(tempNode);
    memset(tempNode, 0, maxsize);
    line.clear();
    

#ifdef PARALLEL
    vector<int>::iterator it;
    int rank;
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    streampos oldpos;
    while (getline(ifs, line)) {
        if (line.compare(0, 1, "*") == 0) {
            ifs.seekg(oldpos);
            break;
        } else {
            size = int(line.size());
            line.copy(tempNode, size, 0);
            tmp = atoi(tempNode) - 1;
            
#ifdef PARALLEL
            if (npart[tmp] == rank) {
                it = find(nodLocal.begin(), nodLocal.end(), tmp);
                if (it == nodLocal.end())
                {
                  ERROR("node %d is not found",tmp);
                }
                else
                {
                  allNodes.push_back(it - nodLocal.begin());
                }
            }
#else
            allNodes.push_back(tmp);        
#endif
            memset(tempNode, 0, maxsize);
        }
        oldpos = ifs.tellg();
    }    
    if (allNodes.size() >0)
    {
        classNodeData* extractD = new classNodeData(name,allNodes,dataType,operation,comp);
        classExtractData::addData(extractD);
    }
};

extern void ElementDataExtractionManagement(vector<classElements *> &elements, ifstream &ifs, vector<int> &epart, vector<int> &elemLocal)
{
    int maxsize = 100;
    char tempNode[maxsize];
    int tmp = 0, size = 0;
    vector<int> commas;
    string line;
    getline(ifs, line);
    tmp = int(line.size());
    int contc = 0;
    for (int j = 0; j < tmp; j++) {
        if (line.compare(j, 1, ",") == 0) {
            commas.push_back(j);
            contc++;
        }
    }
    vector<int> allElemnts;
    bool isAll=false;
    string name, dataType, operation;
    
    
    memset(tempNode, 0, maxsize);
    line.copy(tempNode, commas[0], 0);
    name = tempNode;
    memset(tempNode, 0, maxsize);
    line.copy(tempNode, commas[1] - commas[1 - 1] - 1, commas[1 - 1] + 1);
    dataType = tempNode;
    memset(tempNode, 0, maxsize);
    line.copy(tempNode, line.size(), commas[2 - 1] + 1);
    operation = tempNode;
    memset(tempNode, 0, maxsize);
    line.clear();
    #ifdef PARALLEL
    vector<int>::iterator it;
    int rank;
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    streampos oldpos;
    while (getline(ifs, line)) {
        if (line.compare(0, 1, "*") == 0) {
            ifs.seekg(oldpos);
            break;
        } 
        else if (line.compare(0, 3, "All") == 0){
            isAll = true;
        }
        else {
            size = int(line.size());
            line.copy(tempNode, size, 0);
            tmp = atoi(tempNode) - 1;
            
#ifdef PARALLEL
            if (epart[tmp] == rank) {
                it = find(elemLocal.begin(), elemLocal.end(), tmp);
                if (it == elemLocal.end())
                {
                  ERROR("element %d is not found",tmp);
                }
                else
                {
                  allElemnts.push_back(it - elemLocal.begin());
                }
            }
#else
            allElemnts.push_back(tmp);        
#endif
            memset(tempNode, 0, maxsize);
        }
        oldpos = ifs.tellg();
    }    
    
    if (isAll)
    {
        classElementData* extractD = new classElementData(dataType,operation);
        classExtractData::addData(extractD);
    }
    if (allElemnts.size() >0)
    {
        classElementData* extractD = new classElementData(name,allElemnts,dataType,operation);
        classExtractData::addData(extractD);
    };
};



#ifdef PARALLEL
/*! \brief To be defined (Daniel)
@param[inout] nNod
@param[inout] nElem
@param[inout] NodesPerElement
@param[in] ifs Pointer to the input file (.inp)
*/
extern void infoManagement(int &nNod, int &nElem, int &nodesPerElement, ifstream &ifs) {
    bool flagElement = false;
    bool flagNode = false;
    int tmp2;
    string line, type;
    while (getline(ifs, line)) {
        if (line.compare(0, 5, "*Node") == 0) {
            flagNode = true;
        } else if (line.compare(0, 1, "*") == 0 && flagNode) {
            break;
        } else if (flagNode) {
            nNod++;
        }
    }
    ifs.clear();
    ifs.seekg(0);
    if (!flagNode) {
        ERROR("There is not a list of nodes in the input file");
        exit(-1);
    }
    while (getline(ifs, line)) {
        if (line.compare(0, 8, "*Element") == 0) {
            tmp2 = int(line.size());
            if (line.compare(15, tmp2, "C3D4") == 0) {
                nodesPerElement = 4; // For such type of element
                //         INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "CPS3") == 0) {
                nodesPerElement = 3; // For such type of element
                //          INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "C3D10") == 0) {
                nodesPerElement = 10; // For such type of element
                //          INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "C3D10_lin_surf") == 0) {
                nodesPerElement = 10; // For such type of element
                //          INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "CPS6") == 0) {
                nodesPerElement = 6; // For such type of element
                //          INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "CPS6_lin_surf") == 0) {
                nodesPerElement = 6; // For such type of element
                //          INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "CPE4") == 0) {
                nodesPerElement = 4; // For such type of element
                //          INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "C3D8") == 0) {
                nodesPerElement = 8; // For such type of element
                //          INFO("Type of elements %s",type.c_str());
            } else {
                ERROR("Such type of element is not implemented in MuPhiSim %s", line.c_str());;
                exit(-1);
            }
            flagElement = true;
        } else if (line.compare(0, 1, "*") == 0 && flagElement) {
            break;
        } else if (flagElement) {
            nElem++;
        }
    }
    //Getting the number of elements
    if (!flagElement) {
        ERROR("There is not a list of elements in the input file");
        exit(-1);
    }
}

extern void
infoManagement(int &nNod, int &nElem, int &nodesPerElement, ifstream &ifs, vector<int> &Elems, vector<int> &nodes) {
    bool flagElement = false;
    bool flagNode = false;
    string line, type;
    int maxsize = 100;
    char tmpme2[maxsize];
    memset(tmpme2, 0, maxsize);
    vector<int> commas;
    std::streampos pos;


    while (getline(ifs, line)) {
        if (line.compare(0, 5, "*Node") == 0) {
            flagNode = true;
        } else if (line.compare(0, 1, "*") == 0 && flagNode) {
            break;
        } else if (flagNode) {
            int tmp = int(line.size());
            for (int j = 0; j < tmp; j++) {
                if (line.compare(j, 1, ",") == 0) {
                    commas.push_back(j);
                }
            }
            line.copy(tmpme2, commas[0], 0);
            nodes.push_back(atoi(tmpme2));
            memset(tmpme2, 0, maxsize);
            commas.clear();
            pos = ifs.tellg();
        }
    }
    ifs.clear();
    ifs.seekg(pos);

    if (!flagNode) {
        ERROR("There is not a list of nodes in the input file");
        exit(-1);
    }
    while (getline(ifs, line)) {

        if (line.compare(0, 8, "*Element") == 0) {
            int tmp2 = int(line.size());
            if (line.compare(15, tmp2, "C3D4") == 0) {
                nodesPerElement = 4; // For such type of element
                //         INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "CPS3") == 0) {
                nodesPerElement = 3; // For such type of element
                //          INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "C3D10") == 0) {
                nodesPerElement = 10; // For such type of element
                //          INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "C3D10_lin_surf") == 0) {
                nodesPerElement = 10; // For such type of element
                //          INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "CPS6") == 0) {
                nodesPerElement = 6; // For such type of element
                //          INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "CPS6_lin_surf") == 0) {
                nodesPerElement = 6; // For such type of element
                //          INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "CPE4") == 0) {
                nodesPerElement = 4; // For such type of element
                //          INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "C3D8") == 0) {
                nodesPerElement = 8; // For such type of element
                //          INFO("Type of elements %s",type.c_str());
            }else {
                ERROR("Such type of element is not implemented in MuPhiSim %s", line.c_str());;
                exit(-1);
            }
            flagElement = true;

        } else if (line.compare(0, 1, "*") == 0 && flagElement) {
            break;
        } else if (flagElement) {
            int tmp = int(line.size());
            for (int j = 0; j < tmp; j++) {
                if (line.compare(j, 1, ",") == 0) {
                    commas.push_back(j);
                }
            }
            line.copy(tmpme2, commas[0], 0);
            Elems.push_back(atoi(tmpme2));
            memset(tmpme2, 0, maxsize);
            commas.clear();
            pos = ifs.tellg();

        }
    }
    //Getting the number of elements
    if (!flagElement) {
        ERROR("There is not a list of elements in the input file");
        exit(-1);
    }
    ifs.clear();
    ifs.seekg(pos);
    nNod = nodes.size();
    nElem = Elems.size();

}

extern void FEM_node_elements(int &nNod, int &nElem, ifstream &ifs, vector<int> &elementsFEM, vector<int> &nodesFEM,
                              vector<int> &allElems, vector<int> &allNodes) {
    int maxsize = 100;
    char tempNode[maxsize];
    memset(tempNode, 0, maxsize);
    bool flagNode = false;
    bool flagElem = false;
    string line;
    int tmp = 0, size = 0;

    while (getline(ifs, line)) {

        if (line.compare(0, 13, "*FEM Elements") == 0) {
            flagElem = true;
            flagNode = false;
        } else if (line.compare(0, 10, "*FEM Nodes") == 0) {
            flagNode = true;
            flagElem = false;
        } else {
            if (line.compare(0, 1, "*") == 0 && (flagElem || flagNode)) {
                break;
            } else {
                if (line.compare(0, 3, "All") == 0 && flagElem) {
                    elementsFEM = allElems;
                } else if (line.compare(0, 3, "All") == 0 && flagNode) {
                    nodesFEM = allNodes;
                } else {
                    if (flagNode) {
                        size = int(line.size());
                        line.copy(tempNode, size, 0);
                        tmp = atoi(tempNode) - 1;
                        nodesFEM.push_back(tmp);
                        memset(tempNode, 0, maxsize);
                    } else if (flagElem) {
                        size = int(line.size());
                        line.copy(tempNode, size, 0);
                        tmp = atoi(tempNode) - 1;
                        elementsFEM.push_back(tmp);

                        memset(tempNode, 0, maxsize);
                    }

                }
            }
        }
    }
    nElem = elementsFEM.size();
    nNod = nodesFEM.size();
    if (!flagNode && !flagElem) {
        WARNING("There is not a list of FEM elements or nodes in the input file");
    }
}

extern void FEM_MM_node_elements(ifstream &ifs, vector<int> &elementsFEM, vector<int> &nodesFEM, vector<int> &allElems,
                                 vector<int> &allNodes, vector<int> &nodesMM) {
    int maxsize = 100;
    char tempNode[maxsize];
    memset(tempNode, 0, maxsize);
    bool flagNode = false;
    bool flagElem = false;
    string line;
    int tmp = 0, size = 0;
    bool flagNodeMM = false;
    while (getline(ifs, line)) {

        if (line.compare(0, 13, "*FEM Elements") == 0) {
            flagElem = true;
            flagNode = false;
            flagNodeMM = false;
        } else if (line.compare(0, 10, "*FEM Nodes") == 0) {
            flagNode = true;
            flagElem = false;
            flagNodeMM = false;
        } else if (line.compare(0, 9, "*MM Nodes") == 0) {
            flagNodeMM = true;
            flagNode = false;
            flagElem = false;
        } else {
            if (line.compare(0, 1, "*") == 0 && (flagElem || flagNode || flagNodeMM)) {
                break;
            } else {
                if (line.compare(0, 3, "All") == 0 && flagElem) {
                    elementsFEM = allElems;
                } else if (line.compare(0, 3, "All") == 0 && flagNode) {
                    nodesFEM = allNodes;
                } else if (line.compare(0, 3, "All") == 0 && flagNodeMM) {
                    nodesMM = allNodes;
                } else {
                    if (flagNode) {
                        size = int(line.size());
                        line.copy(tempNode, size, 0);
                        tmp = atoi(tempNode) - 1;
                        nodesFEM.push_back(tmp);
                        memset(tempNode, 0, maxsize);
                    } else if (flagElem) {
                        size = int(line.size());
                        line.copy(tempNode, size, 0);
                        tmp = atoi(tempNode) - 1;
                        elementsFEM.push_back(tmp);

                        memset(tempNode, 0, maxsize);
                    } else if (flagNodeMM) {
                        size = int(line.size());
                        line.copy(tempNode, size, 0);
                        tmp = atoi(tempNode) - 1;
                        nodesMM.push_back(tmp);
                        memset(tempNode, 0, maxsize);
                    }

                }
            }
        }
    }
    if (!flagNode && !flagElem) {
        WARNING("There is not a list of FEM elements or nodes in the input file");
    }
}

/*! \brief To be defined (Daniel)
  @param[inout] elements
  @param[in] ifs Pointer to the input file (.inp)
  @param[inout] nodes
*/
#ifdef METIS
extern void metiselementsManagement(idx_t *elements, ifstream &ifs, idx_t *nodes) {

    int maxsize = 100;
    int index = 0;
    int indexmap = 0;
    bool flagElement = false;
    char tmpme2[maxsize], tempNode[maxsize];
    memset(tempNode, 0, maxsize);
    memset(tmpme2, 0, maxsize);

    string line, type;
    int tmp = 0, tmp2 = 0, numNodes = 0;
    vector<int> myNodes;
    vector<int> commas;
    //classElements *tempElement;
    elements[0] = 0;
    indexmap++;

    //Reading the element information from the file and saving it in temporal variables
    while (getline(ifs, line)) {
        if (line.compare(0, 8, "*Element") == 0) {
            tmp2 = int(line.size());
            if (line.compare(15, tmp2, "C3D4") == 0) {
                type.assign("C3D4");
                numNodes = 4; // For such type of element
                //  INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "CPS3") == 0) {
                type.assign("CPS3");
                numNodes = 3; // For such type of element
                // INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "C3D10") == 0) {
                type.assign("C3D10");
                numNodes = 10; // For such type of element
                // INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "C3D10_lin_surf") == 0) {
                type.assign("C3D10");
                numNodes = 10; // For such type of element
                // INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "CPS6") == 0) {
                type.assign("CPS6");
                numNodes = 6; // For such type of element
                // INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "CPS6_lin_surf") == 0) {
                type.assign("CPS6");
                numNodes = 6; // For such type of element
                // INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "CPE4") == 0) {
                type.assign("CPE4");
                numNodes = 4; // For such type of element
                // INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "C3D8") == 0) {
                type.assign("C3D8");
                numNodes = 8; // For such type of element
                // INFO("Type of elements %s",type.c_str());
            }else {
                ERROR("Such type of element is not implemented in MuPhiSim %s", line.c_str());;
                exit(-1);
            }
            flagElement = true;
        } else {
            if (line.compare(0, 1, "*") == 0 && flagElement) {
                break;
            } else {
                if (flagElement) {
                    tmp = int(line.size());
                    int contc = 0;
                    for (int j = 0; j < tmp; j++) {
                        if (line.compare(j, 1, ",") == 0) {
                            commas.push_back(j);
                            contc++;
                        }
                    }
                    line.copy(tmpme2, commas[0], 0);

                    for (int k = 1; k < numNodes; k++) {
                        memset(tempNode, 0, maxsize);
                        line.copy(tempNode, commas[k] - commas[k - 1] - 1, commas[k - 1] + 1);
                        myNodes.push_back(atoi(tempNode) - 1);
                        memset(tempNode, 0, maxsize);
                    }
                    line.copy(tempNode, tmp - commas[numNodes - 1] - 1, commas[numNodes - 1] + 1);
                    myNodes.push_back(atoi(tempNode) - 1);
                    //Saving the element positioning vector and the nodal structure
                    for (int k = 0; k < numNodes; k++) {
                        nodes[index] = myNodes[k];
                        index++;
                    }
                    elements[indexmap] = elements[indexmap - 1] + numNodes;
                    indexmap++;

                    memset(tmpme2, 0, maxsize);
                    commas.clear();
                    myNodes.clear();
                }
            }
        }
    }

    if (!flagElement) {
        ERROR("There is not a list of elements in the input file");
        exit(-1);
    }
}

//in use
extern void
metiselementsManagement(idx_t *elements, ifstream &ifs, idx_t *nodes, vector<int> elemFEM, vector<int> &mapFEM2global,
                        int nNodeFEM) {

    int maxsize = 100;
    int index = 0;
    int indexmap = 0;
    bool flagElement = false;
    char tmpme2[maxsize], tempNode[maxsize];
    memset(tempNode, 0, maxsize);
    memset(tmpme2, 0, maxsize);

    string line, type;
    int tmp = 0, tmp2 = 0, numNodes = 0;
    vector<int> myNodes;
    vector<int> commas;
    //classElements *tempElement;
    elements[0] = 0;
    indexmap++;
    int countElemNum = 0;
    //Reading the element information from the file and saving it in temporal variables
    while (getline(ifs, line)) {
//        std::cout<<line<<endl;

        if (line.compare(0, 8, "*Element") == 0) {
            tmp2 = int(line.size());
            if (line.compare(15, tmp2, "C3D4") == 0) {
                type.assign("C3D4");
                numNodes = 4; // For such type of element
                //  INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "CPS3") == 0) {
                type.assign("CPS3");
                numNodes = 3; // For such type of element
                // INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "C3D10") == 0) {
                type.assign("C3D10");
                numNodes = 10; // For such type of element
                // INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "C3D10_lin_surf") == 0) {
                type.assign("C3D10");
                numNodes = 10; // For such type of element
                // INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "CPS6") == 0) {
                type.assign("CPS6");
                numNodes = 6; // For such type of element
                // INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "CPS6_lin_surf") == 0) {
                type.assign("CPS6");
                numNodes = 6; // For such type of element
                // INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "CPE4") == 0) {
                type.assign("CPE4");
                numNodes = 4; // For such type of element
                // INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "C3D8") == 0) {
                type.assign("C3D8");
                numNodes = 8; // For such type of element
                // INFO("Type of elements %s",type.c_str());
            }else {
                ERROR("Such type of element is not implemented in MuPhiSim %s", line.c_str());;
                exit(-1);
            }
            flagElement = true;
        } else {
            if (line.compare(0, 1, "*") == 0 && flagElement) {
                break;
            } else {
                if (flagElement) {
                    bool exists = std::find(std::begin(elemFEM), std::end(elemFEM), countElemNum) != std::end(elemFEM);

                    tmp = int(line.size());
                    int contc = 0;
                    for (int j = 0; j < tmp; j++) {
                        if (line.compare(j, 1, ",") == 0) {
                            commas.push_back(j);
                            contc++;
                        }
                    }
                    line.copy(tmpme2, commas[0], 0);

                    for (int k = 1; k < numNodes; k++) {
                        memset(tempNode, 0, maxsize);
                        line.copy(tempNode, commas[k] - commas[k - 1] - 1, commas[k - 1] + 1);
                        myNodes.push_back(atoi(tempNode) - 1);
                        memset(tempNode, 0, maxsize);
                    }
                    line.copy(tempNode, tmp - commas[numNodes - 1] - 1, commas[numNodes - 1] + 1);
                    myNodes.push_back(atoi(tempNode) - 1);


                    if (exists) {
                        //Saving the element positioning vector and the nodal structure
                        for (int k = 0; k < numNodes; k++) {
                            nodes[index] = myNodes[k];
                            index++;
                        }

                        elements[indexmap] = elements[indexmap - 1] + numNodes;
                        mapFEM2global.at(indexmap - 1) = countElemNum;
                        indexmap++;

                    }

                    memset(tmpme2, 0, maxsize);
                    commas.clear();
                    myNodes.clear();
                    countElemNum++;


                }
            }
        }
    }

    bool flagFinish = false;
    int countMatch = 0;

    int k = 0;

    while (!flagFinish) {
        countMatch = 0;

        for (int i = 0; i < nNodeFEM; i++) {
            vector<int> nodesVec(nodes, nodes + elemFEM.size() * numNodes);

            bool nodeExists = std::find(std::begin(nodesVec), std::end(nodesVec), i) != std::end(nodesVec);
            if (!nodeExists) {
                for (int j = 0; j < nodesVec.size(); j++) {
                    if (nodes[j] > i) {
                        nodes[j] = nodes[j] - 1;
                    }
                }
                flagFinish = false;

                countMatch++;
            }

            ++k;
        }
        if (countMatch == 0) {
            flagFinish = true;
            break;
        }
    }
    if (!flagElement) {
        ERROR("There is not a list of elements in the input file");
        exit(-1);
    }

}

extern void metiselementsManagement(idx_t *elements, ifstream &ifs, idx_t *nodes, vector<int> elemFEM, vector<int> &mapFEM2global,
                        int nNodeFEM, vector<int> &mapMM2global, idx_t *elementsMM, idx_t *nodesMM, int nNodeMM) {

    int maxsize = 100;
    int index = 0;
    int indexmap = 0;
    bool flagElement = false;
    char tmpme2[maxsize], tempNode[maxsize];
    memset(tempNode, 0, maxsize);
    memset(tmpme2, 0, maxsize);

    string line, type;
    int tmp = 0, tmp2 = 0, numNodes = 0;
    vector<int> myNodes;
    vector<int> commas;
    //classElements *tempElement;
    elements[0] = 0;
    indexmap++;
    int countElemNum = 0;
    int idxMM_node = 0;
    int idxMM = 1;
    elementsMM[0] = 0;
    //Reading the element information from the file and saving it in temporal variables
    while (getline(ifs, line)) {
        if (line.compare(0, 8, "*Element") == 0) {
            tmp2 = int(line.size());
            if (line.compare(15, tmp2, "C3D4") == 0) {
                type.assign("C3D4");
                numNodes = 4; // For such type of element
                //  INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "CPS3") == 0) {
                type.assign("CPS3");
                numNodes = 3; // For such type of element
                // INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "C3D10") == 0) {
                type.assign("C3D10");
                numNodes = 10; // For such type of element
                // INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "C3D10_lin_surf") == 0) {
                type.assign("C3D10");
                numNodes = 10; // For such type of element
                // INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "CPS6") == 0) {
                type.assign("CPS6");
                numNodes = 6; // For such type of element
                // INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "CPS6_lin_surf") == 0) {
                type.assign("CPS6");
                numNodes = 6; // For such type of element
                // INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "CPE4") == 0) {
                type.assign("CPE4");
                numNodes = 4; // For such type of element
                // INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "C3D8") == 0) {
                type.assign("C3D8");
                numNodes = 8; // For such type of element
                // INFO("Type of elements %s",type.c_str());
            }else {
                ERROR("Such type of element is not implemented in MuPhiSim %s", line.c_str());;
                exit(-1);
            }
            flagElement = true;
        } else {
            if (line.compare(0, 1, "*") == 0 && flagElement) {
                break;
            } else {
                if (flagElement) {
                    bool exists = std::find(std::begin(elemFEM), std::end(elemFEM), countElemNum) != std::end(elemFEM);
                    tmp = int(line.size());
                    int contc = 0;
                    for (int j = 0; j < tmp; j++) {
                        if (line.compare(j, 1, ",") == 0) {
                            commas.push_back(j);
                            contc++;
                        }
                    }
                    line.copy(tmpme2, commas[0], 0);

                    for (int k = 1; k < numNodes; k++) {
                        memset(tempNode, 0, maxsize);
                        line.copy(tempNode, commas[k] - commas[k - 1] - 1, commas[k - 1] + 1);
                        myNodes.push_back(atoi(tempNode) - 1);
                        memset(tempNode, 0, maxsize);
                    }
                    line.copy(tempNode, tmp - commas[numNodes - 1] - 1, commas[numNodes - 1] + 1);
                    myNodes.push_back(atoi(tempNode) - 1);
                    if (exists) {
                        //Saving the element positioning vector and the nodal structure
                        for (int k = 0; k < numNodes; k++) {
                            nodes[index] = myNodes[k];
                            index++;
                        }

                        elements[indexmap] = elements[indexmap - 1] + numNodes;
                        mapFEM2global.at(indexmap - 1) = countElemNum;
                        indexmap++;

                    } else {
                        for (int k = 0; k < numNodes; k++) {
                            nodesMM[idxMM_node] = myNodes[k];
                            idxMM_node++;
                        }

                        elementsMM[idxMM] = elementsMM[idxMM - 1] + numNodes;
                        mapMM2global.at(idxMM - 1) = countElemNum;
                        idxMM++;
                    }

                    memset(tmpme2, 0, maxsize);
                    commas.clear();
                    myNodes.clear();
                    countElemNum++;


                }
            }
        }
    }

    resequence(nodes, nNodeFEM, elemFEM.size() * numNodes);
    resequence(nodesMM, nNodeMM, mapMM2global.size() * numNodes);

    if (!flagElement) {
        ERROR("There is not a list of elements in the input file");
        exit(-1);
    }

}

void resequence(idx_t *nodes, int nodeNum, int vecLength) {
    bool flagFinish = false;
    int countMatch = 0;
    while (!flagFinish) {
        countMatch = 0;

        for (int i = 0; i < nodeNum; i++) {
            vector<int> nodesVec(nodes, nodes + vecLength);

            bool nodeExists = std::find(std::begin(nodesVec), std::end(nodesVec), i) != std::end(nodesVec);
            if (!nodeExists) {

                for (int j = 0; j < nodesVec.size(); j++) {
                    if (nodes[j] > i) {
                        nodes[j] = nodes[j] - 1;
                    }
                }
                flagFinish = false;

                countMatch++;
            }

        }
        if (countMatch == 0) {
            flagFinish = true;
            break;
        }
    }
}

#endif


/*! \brief To be defined (Daniel)
  @param[out] epart Processor ownership mapping of the elements in the whole model epart[element i]=rank (empty for sequential simulations)
  @param[out] npartShare Processor ownership of the nodes in the whole model npartShare[node i]=1 if the nodes is mine or shared with other processor, 0 otherwise (empty for sequential simulations)
  @param[out] nodLocal Mapping of local nodes to global nodes numbering nodlocal[local node index]= global node index (empty for sequential simulations)
  @param[in] ifs Pointer to the input file (.inp)
*/
extern void nodeOwnership(vector<int> &epart, vector<int> &npartShare, vector<int> &nodLocal, ifstream &ifs) {

    int maxsize = 100;
    bool flagElement = false;
    char tmpme2[maxsize], tempNode[maxsize];
    memset(tempNode, 0, maxsize);
    memset(tmpme2, 0, maxsize);
    int rank;
    PetscErrorCode ierr;

    string line, type;
    int tmp = 0, tmp2 = 0, numNodes = 0;
    vector<int> nod_aux;
    vector<int> commas;

    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    CHKERRV(ierr);
    
#ifdef NDEBUG
    // do nothing
#else
    INFO("List of processors assign to each element");
    std::cout << "epart rank" << rank << endl;
    for (int i = 0; i < epart.size(); i++) {
        std::cout << epart[i] << ",";
    }
    std::cout << endl;
#endif 

    //Getting the nodes of my processor and saving it in nod_aux
    while (getline(ifs, line)) {
        if (line.compare(0, 8, "*Element") == 0) {
            tmp2 = int(line.size());
            if (line.compare(15, tmp2, "C3D4") == 0) {
                type.assign("C3D4");
                numNodes = 4; // For such type of element
                //INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "CPS3") == 0) {
                type.assign("CPS3");
                numNodes = 3; // For such type of element
                //	INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "C3D10") == 0) {
                type.assign("C3D10");
                numNodes = 10; // For such type of element
                //	INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "C3D10_lin_surf") == 0) {
                type.assign("C3D10");
                numNodes = 10; // For such type of element
                //	INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "CPS6") == 0) {
                type.assign("CPS6");
                numNodes = 6; // For such type of element
                //	INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "CPS6_lin_surf") == 0) {
                type.assign("CPS6");
                numNodes = 6; // For such type of element
                //	INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "CPE4") == 0) {
                type.assign("CPE4");
                numNodes = 4; // For such type of element
                //	INFO("Type of elements %s",type.c_str());
            } else if (line.compare(15, tmp2, "C3D8") == 0) {
                type.assign("C3D8");
                numNodes = 8; // For such type of element
                //	INFO("Type of elements %s",type.c_str());
            }else {
                ERROR("Such type of element is not implemented in MuPhiSim %s", line.c_str());;
                exit(-1);
            }
            flagElement = true;
        } else {
            if (line.compare(0, 1, "*") == 0 && flagElement) {
                break;
            } else {
                if (flagElement) {

                    tmp = int(line.size());
                    int contc = 0;
                    for (int j = 0; j < tmp; j++) {
                        if (line.compare(j, 1, ",") == 0) {
                            commas.push_back(j);
                            contc++;
                        }
                    }
                    line.copy(tmpme2, commas[0], 0);
                    if (epart[atoi(tmpme2) - 1] == rank) {
                        for (int k = 1; k < numNodes; k++) {
                            memset(tempNode, 0, maxsize);
                            line.copy(tempNode, commas[k] - commas[k - 1] - 1, commas[k - 1] + 1);
                            nod_aux.push_back(atoi(tempNode) - 1);
                            memset(tempNode, 0, maxsize);
                        }
                        line.copy(tempNode, tmp - commas[numNodes - 1] - 1, commas[numNodes - 1] + 1);
                        nod_aux.push_back(atoi(tempNode) - 1);

                    }
                    memset(tmpme2, 0, maxsize);
                    commas.clear();
                }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    //    //Fill with ones the global nodal vector that is own by the processor
    for (int i = 0; (i < nod_aux.size()); i++) {
        npartShare[nod_aux[i]] = 1;
    }

    for (int i = 0; (i < (npartShare.size())); i++) {
        if (npartShare[i] == 1) {
            nodLocal.push_back(i);
        }
    }

}

//Sylvin, make a brief

void readPredef_part(ifstream &ifs, vector<int> &npart, vector<int> &epart) //no label for F or M
{
    bool flagNodes = false;
    bool flagEles = false;

    fill(npart.begin(), npart.end(), 99);
    fill(epart.begin(), epart.end(), 99);

    string line, nextLine;
    int maxsize = 100;
    char tempChar[maxsize];
    memset(tempChar, 0, maxsize);
    int tmp = 0, size = 0, proc = 0;

    while (getline(ifs, line)) {
        if (line.compare(0, 6, "*ELSET") == 0) {
            flagEles = true;
            flagNodes = false;

            getline(ifs, nextLine);
            if (nextLine.compare(0, 5, "*CPU=") == 0) {
                size = int(nextLine.size());
                nextLine.copy(tempChar, size - 5, 5);
                proc = atoi(tempChar) - 1;
            } else {
                ERROR("Processor number is not assigned for Elset");
                exit(-1);
            }

            memset(tempChar, 0, maxsize);
        } else if (line.compare(0, 5, "*NSET") == 0) {
            flagNodes = true;
            flagEles = false;
            getline(ifs, nextLine);
//            std::cout<<nextLine<<endl;
            if (nextLine.compare(0, 5, "*CPU=") == 0) {
                size = int(nextLine.size());
                nextLine.copy(tempChar, size - 5, 5);
                proc = atoi(tempChar) - 1;
//                std::cout<<"CPU get "<<proc<<endl;
            } else {
                ERROR("Processor number is not assigned for Nset");
                exit(-1);
            }
            memset(tempChar, 0, maxsize);
        } else {
            if (line.compare(0, 1, "*") == 0 && (flagEles || flagNodes)) {
//                std::cout<<"end"<<endl;
                break;
            } else {
                if (flagNodes) {
                    size = int(line.size());
                    line.copy(tempChar, size, 0);
                    tmp = atoi(tempChar) - 1;
                    npart[tmp] = proc;
                    memset(tempChar, 0, maxsize);
                } else if (flagEles) {
                    size = int(line.size());
                    line.copy(tempChar, size, 0);
                    tmp = atoi(tempChar) - 1;
                    epart[tmp] = proc;
                    memset(tempChar, 0, maxsize);
                }
            }
        }
    }
}

//Sylvin, make a brief

bool meshPartitionMethod(ifstream &ifs) {
    string line;
    bool flagMETIS = true;
    ifs.seekg(0);
    while (getline(ifs, line)) {
        if ((line.compare(0, 6, "*ELSET") == 0) || (line.compare(0, 5, "*NSET") == 0)) {
            flagMETIS = false;
            break;
        }
    }
    return flagMETIS;
}

#endif
