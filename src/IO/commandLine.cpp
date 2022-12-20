//
//
// File authors: see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file commandLine.cpp
  \brief Functions related to the management of the flags given in the command line.
*/
#include "commandLine.h"
#include <sys/stat.h>
#ifdef METIS
#include <metis.h>
#endif

std::string GeneralOptions::workingDirectoryName="";
std::string GeneralOptions::inputDirectoryName="input"; // defaut location
std::string GeneralOptions::outputDirectoryName="output"; // defaut location
std::string GeneralOptions::inputFileName="";
std::string GeneralOptions::petscOptionsFile="options.petsc";
int GeneralOptions::commRank = 0;
int GeneralOptions::commSize = 1;
double GeneralOptions::MMScalingFactor = 1.;
int GeneralOptions::MMIntegrationOrder=0; // defaut intgration order of background mesh, 
                                        // a value = 0 means defaut by each element, 
                                        // a value >0 will applied for all elements
bool GeneralOptions::noImplicitIterative=false; // true if the implicit is performed iteratively,  

bool GeneralOptions::MMAdaptiveSupportRadius=true; // true if the support radius is increased by a factor 1.1 
                                                    // if the optimisation problem to find the shape functions 
                                                    // cannot be successfully performed, see shapeFunctions.h 

void printToScreen()
{
    printf("MuPhiSim: a platform for multiphysics simulations\n");
    printf("Usage: MuPhiSim-executable [options] [value]\n");
    printf("Related options name if exist, are given below:\n");
    printf("-input string\t Set the name of input file\n");
    printf("-inputDir string\t Set the input folder containing the input file [by default: input]\n");
    printf("-outputDir string\t Set the output folder [by default: output]\n");
    printf("-MMScalingFactor value\t Set the meshless scaling factor by a double value [by default: 1]\n");
    printf("-MMIntegOrder n\t Set the integration order of the meshless background mesh by a positive integer (0,1,2,...),\n\t [by default: 0; the order decided by element type]\n");
    printf("-MMAdaptiveSupportRadius n\t Activate (1 is used)/desativate (0 is used) the adaptive support radius in meshless simulation [by default: 1]\n");
    printf("-noImplicitIterative n\t Activate (1 is used)/desativate (0 is used) the iterative process in the implicit solver [by default: 0]\n");
    printf("-petscOptionsFile fileName \t Using options for petsc with file [by default options.petsc in the current folder]\n");
};

void GeneralOptions::initialise(int argc, char **argv)
{
    char cwd[1024];
    string str("");
    if (getcwd(cwd, sizeof(cwd)) == NULL) {
        ERROR("The program cannot define the working directory");
        exit(-1);
    };
    workingDirectoryName.append(cwd);
    if (argc==1)
    {
        printToScreen();
        exit(0);
    }
    for (int i = 1; i < argc; i++) 
    {
        if ((strcmp(argv[i], "-h")==0) || (strcmp(argv[i], "-help")==0))
        {
            printToScreen();
            exit(0);
        }
        else if (strcmp(argv[i], "-input")==0) 
        {
            inputFileName="";
            inputFileName.append(argv[i+1]);
        }
        else if (strcmp(argv[i], "-inputDir")==0)
        {
            inputDirectoryName = "";
            inputDirectoryName.append(argv[i+1]);
        }
        else if (strcmp(argv[i], "-outputDir")==0)
        {
            outputDirectoryName = "";
            outputDirectoryName.append(argv[i+1]);
        }
        else if (strcmp(argv[i], "-MMScalingFactor")==0)
        {
            GeneralOptions::MMScalingFactor = atof(argv[i+1]);
            INFO("Meshless scaling factor: %g",GeneralOptions::MMScalingFactor);            
        }
        else if (strcmp(argv[i], "-MMIntegOrder")==0)
        {
            GeneralOptions::MMIntegrationOrder = atoi(argv[i+1]);
            INFO("Meshless background mesh integration order: %d",GeneralOptions::MMIntegrationOrder);            
        }
        else if (strcmp(argv[i], "-MMAdaptiveSupportRadius")==0)
        {
            GeneralOptions::MMAdaptiveSupportRadius = (bool)atoi(argv[i+1]);
            if (GeneralOptions::MMAdaptiveSupportRadius)
            {
                INFO("Adaptive meshless support radius is considered");   
            }
            else
            {
                INFO("Adaptive meshless support radius is not considered");  
            }
        }
        else if (strcmp(argv[i], "-noImplicitIterative")==0)
        {
            GeneralOptions::noImplicitIterative = (bool)atoi(argv[i+1]);
            if (GeneralOptions::noImplicitIterative)
            {
                INFO("Implicit scheme is performed one iteration in each time step");
            }
        }
        else if (strcmp(argv[i], "-petscOptionsFile")==0)
        {
            GeneralOptions::petscOptionsFile="";
            GeneralOptions::petscOptionsFile.append(argv[i+1]);
            INFO("petsc options are considered in %s",GeneralOptions::petscOptionsFile.c_str());
        }
    }
    
    // petsc
    struct stat buffer;   
    if (stat (GeneralOptions::petscOptionsFile.c_str(), &buffer) == 0)
    {
        //INFO("petsc options file %s exists",GeneralOptions::petscOptionsFile.c_str());
        PetscInitialize(PETSC_NULL, PETSC_NULL, GeneralOptions::petscOptionsFile.c_str(), PETSC_NULL);
    }
    else
    {
        //INFO("petsc options file %s does not exist",GeneralOptions::petscOptionsFile.c_str());
        PetscInitialize(PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);
    }
    
    
#ifdef PARALLEL
    MPI_Comm_rank(PETSC_COMM_WORLD, &commRank);
    MPI_Comm_size(PETSC_COMM_WORLD, &commSize);
#endif // PARALLEL
    
    if (commRank == 0)
    {
        INFO("Current working directory: %s",workingDirectoryName.c_str());
        INFO("Input file: %s",inputFileName.c_str());
        INFO("Input directory: %s",inputDirectoryName.c_str());
        INFO("Output directory: %s",outputDirectoryName.c_str());
    }    

};

void GeneralOptions::finalise()
{
    PetscFinalize();
};

string GeneralOptions::getInputFilePath()
{
    string str(workingDirectoryName);
    str.append("/"+inputDirectoryName+"/"+inputFileName);
    return str;
};

string GeneralOptions::getOutputDirectoryPath()
{
    string str(workingDirectoryName);
    str.append("/"+outputDirectoryName);
    return str;
};


/*! \brief It checks some logical numbers obtained from the input file.
  @param[in] argc Command line size
  @param[in] argv Command line data
  @param[out] nod Array with all nodes in the domain
  @param[out] elem Array with all elements (or background cells) in the domain
  @param[out] NeumannBCs Array with the Neumann boundary conditions
  @param[out] ndim Dimension of the problem
  @param[out] nodFEM Array with the labels that correspond to the FEM nodes in the nodes general array
  @param[out] nodMM Array with the labels that correspond to the MM nodes in the nodes general array
  @param[out] nodShare Array with the labels that correspond to the shared nodes in the nodes general array
  @param[out] elemFEM Array with the labels that correspond to the FEM elements in the elements general array
  @param[out] elemMM Array with the labels that correspond to the MM background cells in the elements general array
  @param[out] geoShape This array contains all elements to define the geometrical boundary of the domain. Therefore, the elements contained are for a dimension lower than ndim (i.e., line elements for a simulation in 2D and surface elements for 3D)
  @param[out] consMod Array with the constitutive models in the domain
  @param[out] epart Processor ownership mapping of the elements in the whole model epart[element i]=rank (empty for sequential simulations)
  @param[out] npartShare Processor ownership of the nodes in the whole model npartShare[node i]=1 if the nodes is mine or shared with other processor, 0 otherwise (empty for sequential simulations)
  @param[out] nodLocal Mapping of local nodes to global nodes numbering nodlocal[local node index]= global node index (empty for sequential simulations)
  @param[out] elemLocal Mapping of local elements to global elements numbering -> elemlocal[local element index] = global element index (empty for sequential simulations)
  @param[out] npart Processor ownership mapping of the nodes in the whole model npart[node i]=rank (empty for sequential simulations)
  @param[out] solver Array with all solvers to be executed
*/
extern void inputManagement(vector<classNodes *> &nod, vector<classElements *> &elem,
                            vector<classNeumannBCs *> &NeumannBCs, int &ndim, vector<int> &nodFEM, vector<int> &nodMM,
                            vector<int> &nodShare, vector<int> &elemFEM, vector<int> &elemMM,
                            vector<classElements *> &geoShape, vector<constitutiveModels *> &consMod,
                            vector<int> &epart, vector<int> &npartShare, vector<int> &nodLocal, vector<int> &elemLocal,
                            vector<int> &npart, vector<solvers *> &solver, vector<double> &iniValueExtraDof) {
    ifstream ifs;
    int surface_order = -1;
    string str = GeneralOptions::getInputFilePath();
    ifs.open(str.c_str(), ifstream::in);
    if (ifs.is_open()) {
        nodesManagement(nod, ndim, ifs, npartShare);
        if (nod.size() == 0) {
            ERROR("There should be a general list of nodes in the input file");
            exit(-1);
        }
        ifs.clear();
        ifs.seekg(0);
        elementsManagement(elem, ndim, ifs, nod, epart, nodLocal, surface_order);
        if (elem.size() == 0) {
            ERROR("There should be a general list of elements in the input file");
            exit(-1);
        }
        ifs.clear();
        ifs.seekg(0);
        nodesFEMManagement(nodFEM, ifs, nod, npartShare);
        ifs.clear();
        ifs.seekg(0);
        nodesMMManagement(nodMM, ifs, nod, npartShare);
        ifs.clear();
        ifs.seekg(0);
        nodesShareManagement(nodShare, ifs, npartShare);
        ifs.clear();
        ifs.seekg(0);
        int size1 = nodFEM.size();
        int size2 = nodMM.size();
        int size = nodShare.size();
        if (size1 == 0 && size2 == 0) {
            ERROR("The list of nod in the input file should be either FEM or MM");
            exit(-1);
        }
        if (size == 0) {
            if (size1 == 0 && size2 != 0) {
                WARNING("The list of nodes in the input file is assigned to MM");
            }
            if (size1 != 0 && size2 == 0) {
                WARNING("The list of nodes in the input file is assigned to FEM");
            }
        } else {
            if (size1 == 0 || size2 == 0) {
                ERROR("There are shared nodes, but only either FEM or MM have been defined in the input file");
            }
        }
        elementsFEMManagement(elemFEM, ifs, elem, epart);

        ifs.clear();
        ifs.seekg(0);
        elementsMMManagement(elemMM, ifs, elem, epart);

        ifs.clear();
        ifs.seekg(0);
        size1 = elemFEM.size();
        size2 = elemMM.size();
        if (size1 == 0 && size2 == 0) {
            ERROR("There should be either FEM elements or MM elements");
            exit(-1);
        }
        materialManagement(elem, ifs, ndim, consMod, epart, elemLocal);
        if (consMod.size() == 0) {
            ERROR("The constitutive model has not been defined in the input file");
            exit(-1);
        }
        ifs.clear();
        ifs.seekg(0);

        initialConditionsManagement(ifs,
                                    iniValueExtraDof); //(Optional) This function allows to define an initial vale different from zero for extradof, e.g., Temperature
        ifs.clear();
        ifs.seekg(0);

        int nDof = nod.size() * ndim; // Only considering the displacements so far
        solverManagement(solver, ifs, elem, nod, nDof, ndim, npart, nodLocal, epart, elemLocal);

        if (solver.size() == 0) {
            ERROR("The solver has not been defined in the input file");
            exit(-1);
        }
        
        ifs.clear();
        ifs.seekg(0);
        dataExtractionManagement(ifs, elem, nod, npart, nodLocal, epart, elemLocal);
        
        ifs.clear();
        ifs.seekg(0);
        NeumannCommonManagement(elem, nod, NeumannBCs, ndim, ifs, npart, nodLocal, epart, elemLocal,
                                surface_order);

        ifs.clear();
        ifs.seekg(0);
    } else {
        ERROR("The input file cannot be opened");
        exit(-1);
    }

}

//with flag for MM
extern void inputManagement(vector<classNodes *> &nod, vector<classElements *> &elem,
                            vector<classNeumannBCs *> &NeumannBCs, int &ndim, vector<int> &nodFEM, vector<int> &nodMM,
                            vector<int> &nodShare, vector<int> &elemFEM, vector<int> &elemMM,
                            vector<classElements *> &geoShape, vector<constitutiveModels *> &consMod,
                            vector<int> &epart, vector<int> &npartShare, vector<int> &nodLocal, vector<int> &elemLocal,
                            vector<int> &npart, vector<solvers *> &solver, bool flagProc4MM,
                            vector<double> &iniValueExtraDof) {
    ifstream ifs;
    int surface_order = -1;
    string str = GeneralOptions::getInputFilePath();
    ifs.open(str.c_str(), ifstream::in);
    if (ifs.is_open()) {
        nodesManagement(nod, ndim, ifs, npartShare);

        if (nod.size() == 0) {
            ERROR("There should be a general list of nodes in the input file");
            exit(-1);
        }
        ifs.clear();
        ifs.seekg(0);

        elementsManagement(elem, ndim, ifs, nod, epart, nodLocal, surface_order);

        if (elem.size() == 0) {
            ERROR("There should be a general list of elements in the input file");
            exit(-1);
        }
        ifs.clear();
        ifs.seekg(0);
        nodesFEMManagement(nodFEM, ifs, nod, npartShare, nodLocal);


        ifs.clear();
        ifs.seekg(0);
        nodesMMManagement(nodMM, ifs, nod, npartShare, nodLocal);

        ifs.clear();
        ifs.seekg(0);
        nodesShareManagement(nodShare, ifs, npartShare);

        ifs.clear();
        ifs.seekg(0);
        int size1 = nodFEM.size();
        int size2 = nodMM.size();
        int size = nodShare.size();
        if (size1 == 0 && size2 == 0) {
            ERROR("The list of nod in the input file should be either FEM or MM");
            exit(-1);
        }
        if (size == 0) {
            if (size1 == 0 && size2 != 0) {
                WARNING("The list of nodes in the input file is assigned to MM");
            }
            if (size1 != 0 && size2 == 0) {
                WARNING("The list of nodes in the input file is assigned to FEM");
            }
        } else {
            if (size1 == 0 || size2 == 0) {
                ERROR("There are shared nodes, but only either FEM or MM have been defined in the input file");
            }
        }
        elementsFEMManagement(elemFEM, ifs, elem, epart, elemLocal);

        ifs.clear();
        ifs.seekg(0);
        elementsMMManagement(elemMM, ifs, elem, epart, elemLocal);

        ifs.clear();
        ifs.seekg(0);
        size1 = elemFEM.size();
        size2 = elemMM.size();
        if (size1 == 0 && size2 == 0) {
            ERROR("There should be either FEM elements or MM elements");
            exit(-1);
        }
        materialManagement(elem, ifs, ndim, consMod, epart, elemLocal);


        if (consMod.size() == 0) {
            ERROR("The constitutive model has not been defined in the input file");
            exit(-1);
        }
        ifs.clear();
        ifs.seekg(0);

        initialConditionsManagement(ifs,
                                    iniValueExtraDof); //(Optional) This function allows to define an initial vale different from zero for extradof, e.g., Temperature
        ifs.clear();
        ifs.seekg(0);

        int nDof = nod.size() * ndim; // Only considering the displacements so far
        solverManagement(solver, ifs, elem, nod, nDof, ndim, npart, nodLocal, epart, elemLocal);

        if (solver.size() == 0) {
            ERROR("The solver has not been defined in the input file");
            exit(-1);
        }
        
        ifs.clear();
        ifs.seekg(0);
        dataExtractionManagement(ifs, elem, nod, npart, nodLocal, epart, elemLocal);
        
        ifs.clear();
        ifs.seekg(0);
        NeumannCommonManagement(elem, nod, NeumannBCs, ndim, ifs, npart, nodLocal, epart, elemLocal,
                                surface_order);

        ifs.clear();
        ifs.seekg(0);
    } else {
        ERROR("The input file cannot be opened");
        exit(-1);
    }
}

#ifdef PARALLEL
/*! \brief It reads the command line and calls other input related functions (DANIEL)
  @param[in] argc Number of arguments in the command line
  @param[in] argv Array with the command line strings
  @param[out] nod Array with all nodes in the domain
  @param[out] elem Array with all elements (or background cells) in the domain
  @param[out] DirichletNod Array with the Dirichlet boundary nodes
  @param[out] ndim Dimension of the problem
  @param[out] nodFEM Array with the labels that correspond to the FEM nodes in the nodes general array
  @param[out] nodMM Array with the labels that correspond to the MM nodes in the nodes general array
  @param[out] nodShare Array with the labels that correspond to the shared nodes in the nodes general array
  @param[out] elemFEM Array with the labels that correspond to the FEM elements in the elements general array
  @param[out] elemMM Array with the labels that correspond to the MM background cells in the elements general array
  @param[out] scaleFactor This factor will be multiplied by the time step in order to be safer (further away from the critical time step)
  @param[out] geoShape This array will contain all elements to define the geometrical boundary of the domain. Therefore, the elements contained will be defined for a dimension lower than ndim (i.e., line elements for a simulation in 2D and surface elements for 3D)
*/

// This was the previous modelInfo version (June 2019)

// extern void modelInfo(int argc, char **argv, int &nNod, int &nElem, int &nodesPerElement, vector<int> &elementsFEM, vector<int> &nodesFEM, vector<int> &nodesMM){

//     ifstream ifs;
//     char cwd[1024];
//     string str;
//     vector<int> allElems;
//     vector<int> allNodes;
//     if(getcwd(cwd, sizeof(cwd))==NULL){
//         ERROR("The program cannot define the working directory");
//         exit(-1);
//     };
//     str.append(cwd);
//     str.append("/input/");
//     for(int i=1;i<argc;i++){
//         if(strcmp(argv[i],"-input")){
//             str.append(argv[i]);
//             ifs.open(str.c_str(), ifstream::in);
//             if(ifs.is_open()){
//                 infoManagement(nNod, nElem, nodesPerElement, ifs,allElems,allNodes);

//                 FEM_MM_node_elements(ifs,elementsFEM,nodesFEM,allElems,allNodes,nodesMM);
//             }
//             else{
//                 ERROR("The input file cannot be opened");
//                 exit(-1);
//             }
//         }
//     }
// }

// This is the new Dongli's function implemented for partitioning meshless

extern void modelInfo(int &nNod, int &nElem, int &nodesPerElement, vector<int> &elementsFEM,
                      vector<int> &nodesFEM, vector<int> &nodesMM, bool &flagMETIS) {

    ifstream ifs;
    string str = GeneralOptions::getInputFilePath();
    vector<int> allElems;
    vector<int> allNodes;
    ifs.open(str.c_str(), ifstream::in);
    if (ifs.is_open()) {
        infoManagement(nNod, nElem, nodesPerElement, ifs, allElems, allNodes);
        //           string line;
        //           getline(ifs,line);
        //           std::cout<<"model info nNode total\t"<<nNod<<"\n";
        //           FEM_node_elements(nNodeFEM,nElemFEM, ifs,elementsFEM,nodesFEM,allElems,allNodes);
//                std::cout<<"infoManagement"<<endl;
        FEM_MM_node_elements(ifs, elementsFEM, nodesFEM, allElems, allNodes, nodesMM);

        //           std::cout<<"nNodeFEM"<<nNodeFEM<<"\t nElemFEM"<<nElemFEM<<"\t nNode"<<nNod<<endl;
        //           ifs.clear();
        //           ifs.seekg(0);
//                std::cout<<"FEM_MM"<<endl;
        flagMETIS = meshPartitionMethod(ifs);
//                std::cout<<"flagMETIS" <<flagMETIS<<endl;

    } else {
        ERROR("The input file cannot be opened");
        exit(-1);
    }
}

// This is the second mew modelInfo function for FEM Mesh only (Dongli)

extern void modelInfo(int &nNod, int &nElem, int &nodesPerElement, int &nNodeFEM, int &nElemFEM,
                      vector<int> &elementsFEM, vector<int> &nodesFEM) {

    ifstream ifs;
    string str = GeneralOptions::getInputFilePath();
    vector<int> allElems;
    vector<int> allNodes;
    ifs.open(str.c_str(), ifstream::in);
    if (ifs.is_open()) {
        infoManagement(nNod, nElem, nodesPerElement, ifs, allElems, allNodes);
        //           string line;
        //           getline(ifs,line);
        //           std::cout<<"model info nNode total\t"<<nNod<<"\n";
        FEM_node_elements(nNodeFEM, nElemFEM, ifs, elementsFEM, nodesFEM, allElems, allNodes);
        //           std::cout<<"nNodeFEM"<<nNodeFEM<<"\t nElemFEM"<<nElemFEM<<"\t nNode"<<nNod<<endl;
        //           ifs.clear();
        //           ifs.seekg(0);
    } else {
        ERROR("The input file cannot be opened");
        exit(-1);
    }
}
/*! \brief rearrange elements and nodes for metis
  @param[in] argc Command line size
  @param[in] argv Command line data
  @param[inout] nod All associated node numbers for FEM elements
  @param[inout] elem Index of starting node index in the nod array
  @param[in] flagCombined Flag to indicate the FEM/MM mode
  @param[in] ElemFEM Element numbers for FEM elements
  @param[inout] nodeMM Associated node numbers for MM elements
  @param[inout] mapFEM2global Mapping of the indices in FEM to the global arrays
*/

// Dongli's new polymorphism for metis input
#ifdef METIS
extern void metisInput(idx_t *nod, idx_t *elem, bool flagCombined, vector<int> ElemFEM,
                       vector<int> &mapFEM2global, int nNodeFEM) {

    ifstream ifs;
    string str = GeneralOptions::getInputFilePath();
    ifs.open(str.c_str(), ifstream::in);
    if (ifs.is_open()) {
        if (flagCombined) {
            metiselementsManagement(elem, ifs, nod, ElemFEM, mapFEM2global, nNodeFEM);
        } else {
            metiselementsManagement(elem, ifs, nod);
        }
        ifs.clear();
        ifs.seekg(0);
    } else {
        ERROR("The input file cannot be opened");
        exit(-1);
    }

}

// This is the previous with MM parameters
extern void metisInput(idx_t *nod, idx_t *elem, bool flagCombined, vector<int> ElemFEM,
                       vector<int> &mapFEM2global, int nNodeFEM, vector<int> &mapMM2global, idx_t *elementsMM,
                       idx_t *nodesMM, int nNodeMM) {

    ifstream ifs;
    string str = GeneralOptions::getInputFilePath();

    ifs.open(str.c_str(), ifstream::in);
    if (ifs.is_open()) {
        if (flagCombined) {
            metiselementsManagement(elem, ifs, nod, ElemFEM, mapFEM2global, nNodeFEM, mapMM2global, elementsMM,
                                    nodesMM, nNodeMM);
        } else {
            metiselementsManagement(elem, ifs, nod);
        }
        ifs.clear();
        ifs.seekg(0);
    } else {
        ERROR("The input file cannot be opened");
        exit(-1);
    }
}
#endif //METIS

/*! \brief It checks some logical numbers obtained from the input file (DANIEL)
  @param[out] npart Processor ownership mapping of the nodes in the whole model npart[node i]=rank (empty for sequential simulations)
  @param[out] epart Processor ownership mapping of the elements in the whole model epart[element i]=rank (empty for sequential simulations)
  @param[out] nodLocal Local mapping of local nodes to global nodes numbering nodlocal[local node index]= global node index (empty for sequential simulations)
  @param[out] nodLocal Mapping of local nodes to global nodes numbering nodlocal[local node index]= global node index (empty for sequential simulations)
  @param[out] elemLocal Mapping of local elements to global elements numbering -> elemlocal[local element index] = global element index (empty for sequential simulations)
  @param[out] npartShare Processor ownership of the nodes in the whole model npartShare[node i]=1 if the nodes is mine or shared with other processor, 0 otherwise (empty for sequential simulations)
*/
extern void getownership(vector<int> &npart, vector<int> &epart, vector<int> &nodLocal,
                         vector<int> &elemLocal, vector<int> &npartShare) {
    int rank = GeneralOptions::commRank;
    ifstream ifs;
    string str = GeneralOptions::getInputFilePath();
    //Building global to local element map matrix:
    for (int i = 0; i < (epart.size()); i++) {
        if (epart[i] == rank) {
            elemLocal.push_back(i);
        }
    }

    //Building the nodes ownership vector and local nodes map matrix

    ifs.open(str.c_str(), ifstream::in);
    if (ifs.is_open()) {
        nodeOwnership(epart, npartShare, nodLocal, ifs);
        ifs.clear();
        ifs.seekg(0);
    } else {
        ERROR("The input file cannot be opened");
        exit(-1);
    }
}


/*! \brief It checks some logical numbers obtained from the input file.
  @param[out] nNod Number of global nodes
  @param[out] nNod Number of nodes
  @param[out] nElem number of Elements
  @param[out] epart Processor ownership mapping of the elements in the whole model epart[element i]=rank (empty for sequential simulations)
  @param[out] npartShare Processor ownership of the nodes in the whole model npartShare[node i]=1 if the nodes is mine or shared with other processor, 0 otherwise (empty for sequential simulations)
  @param[out] nodLocal Mapping of local nodes to global nodes numbering nodlocal[local node index]= global node index (empty for sequential simulations)
  @param[out] elemLocal Mapping of local elements to global elements numbering -> elemlocal[local element index] = global element index (empty for sequential simulations)
  @param[out] npart Processor ownership mapping of the nodes in the whole model npart[node i]=rank (empty for sequential simulations)
*/

// This was the previous parallelManagement function (June 2019)

// extern void parallelManagement(int argc, char **argv, int &nNodGlobal, int &nNod, int &nElem, vector<int> &epart, vector<int> &npartShare,vector<int> &nodLocal,vector<int> &elemLocal, vector<int> &npart,bool &flagProc4MM)
// {
//     int nprcs=1; // MPI number of processor and total number of processors, respectively
//     int nodesPerElement;
//     nprcs=0;
//     int rank=0;
//     //Recover rank and number of processes information
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);//CHKERRQ(ierr);
//     MPI_Comm_size(PETSC_COMM_WORLD,&nprcs);//CHKERRQ(ierr);
//     int nNodeFEM = 0;
//     int nElemFEM = 0;
//     int nNodeMM = 0;
//     vector<int> elementsFEM;
//     vector<int> nodesFEM;
//     vector<int> nodesMM;

//     bool flagCombined = false;
//     int procNumMM = 1;

//     if (rank==0)
//     {
//         INFO("Running simulation on %i cores", nprcs);
//         INFO("Partitioning mesh ...");
//         //Obtain the model information from the input file (!!! In the future just need inputMetis(elements[nelem,nnodes_perelement])
//         modelInfo(argc,argv,nNod,nElem,nodesPerElement,elementsFEM,nodesFEM,nodesMM);
//         nNodGlobal=nNod;
//         nNodeFEM = nodesFEM.size();
//         nElemFEM = elementsFEM.size();
//         nNodeMM = nodesMM.size();
//         int nElemMM = nElem - nElemFEM;
//         int GPsPerMM = 6;
//         int GPsPerFEM = 1;// XXX need to be output directly from the element type
//         procNumMM = int(nprcs * nElemMM * GPsPerMM/(nElemMM * GPsPerMM+nElemFEM*GPsPerFEM));//XXX
//     }

//     //Casting the information from the root process
//     MPI_Barrier(MPI_COMM_WORLD); //Wait all the processors to reach this point
//     MPI_Bcast(&nNod,1,MPI_INT,0,MPI_COMM_WORLD);
//     MPI_Bcast(&nNodGlobal,1,MPI_INT,0,MPI_COMM_WORLD);
//     MPI_Bcast(&nElem,1,MPI_INT,0,MPI_COMM_WORLD);
//     MPI_Bcast(&nodesPerElement,1,MPI_INT,0,MPI_COMM_WORLD);
//     MPI_Bcast(&nNodeFEM,1,MPI_INT,0,MPI_COMM_WORLD);
//     MPI_Bcast(&nElemFEM,1,MPI_INT,0,MPI_COMM_WORLD);

//     //Creating the METIS variables
//     idx_t nElemMetis=nElemFEM;//nElem;                               // Number of elements in the global model
//     idx_t nNodMetis=nNodeFEM; //nNod;                                 // Number of nodes in the global model
//     idx_t nElemMetisMM=nElem - nElemFEM;
//     idx_t nNodMetisMM=nNodeMM;
//     idx_t  nprcsMetis=nprcs-procNumMM;                              // default for combined FEM MM
//     idx_t  nprcsMetisMM=procNumMM;


//     vector<int> mapMM2global(nElem - nElemFEM,0);

//     int maxProcNum = 0;//nprcsMetis-procNumMM;

//     if (nNodeFEM ==nNod)
//     {
//         nprcsMetis=nprcs;                              // FEM only
//     }
//     else if(nNodeFEM == 0) //MM only
//     {
//         // for MM

//         nElemMetis=nElem - nElemFEM;
//         nNodMetis=nNodeMM;
//         nprcsMetis=nprcs;                              // FEM only
//         flagProc4MM = true;
//     }
//     else
//     {
//         flagCombined = true;
//         nElemMetisMM=nElem - nElemFEM;
//         nNodMetisMM=nNodeMM;
//     }

//     // Partitioning of the mesh by Metis
//     idx_t npartMetis[nNodMetis];                         // Processor mapping of the nodes
//     assert(npartMetis != NULL);

//     idx_t  nMetis;                                        // Number of Parts in Metis
//     idx_t  epartMetis[nElemMetis];                       // Processor mappingof th elements
//     assert(epartMetis != NULL);
//     vector<int> mapFEM2global(nElemFEM,0);
//     idx_t nparttmp[nNod];
//     idx_t eparttmp[nElem];

//     // for MM
//     MPI_Bcast(&nNodeMM,1,MPI_INT,0,MPI_COMM_WORLD);

//     idx_t npartMetisMM[nNodMetisMM];                         // Processor mapping of the nodes
//     assert(npartMetisMM != NULL);

//     idx_t  nMetisMM;                                        // Number of Parts in Metis
//     idx_t  epartMetisMM[nElemMetisMM];                       // Processor mappingof th elements
//     assert(epartMetisMM != NULL);

//     vector<int>::iterator it;

//     if (rank==0)
//     {
//         idx_t elemetisMetis[nodesPerElement*nElemMetis];  // Vector of nodes of all the elements
//         assert(elemetisMetis != NULL);
//         idx_t mapmetisMetis[nElemMetis+1];  // Pre-metis vector distribution of elements in elemetisMetis
//         assert(elemetisMetis != NULL);
//         //MM
//         idx_t elemetisMetisMM[nodesPerElement*nElemMetisMM];  // Vector of nodes of all the elements
//         assert(elemetisMetisMM != NULL);
//         idx_t mapmetisMetisMM[nElemMetisMM+1];  // Pre-metis vector distribution of elements in elemetisMetis
//         assert(elemetisMetisMM != NULL);

//         metisInput(argc,argv,elemetisMetis,mapmetisMetis, flagCombined, elementsFEM,mapFEM2global,nNodeFEM,mapMM2global,mapmetisMetisMM,elemetisMetisMM,nNodeMM);

//         if(flagCombined)
//         {
//             METIS_PartMeshNodal(&nElemMetis,&nNodMetis,mapmetisMetis,elemetisMetis,NULL,NULL,&nprcsMetis,NULL,NULL,&nMetis,epartMetis,npartMetis);
//             METIS_PartMeshNodal(&nElemMetisMM,&nNodMetisMM,mapmetisMetisMM,elemetisMetisMM,NULL,NULL,&nprcsMetisMM,NULL,NULL,&nMetisMM,epartMetisMM,npartMetisMM);

//             vector<int> tmpConv(epartMetis, epartMetis + nElemMetis);
//             for(int k=0; k<tmpConv.size(); k++)
//             {
//                 if(tmpConv[k]>maxProcNum)
//                 {
//                     maxProcNum = tmpConv[k];
//                 }
//             }

//             for(int i =0; i<nElem;i++)
//             {
//                 bool exists = std::find(std::begin(mapFEM2global), std::end(mapFEM2global), i) != std::end(mapFEM2global);

//                 if(exists)
//                 {
//                     it=find(mapFEM2global.begin(), mapFEM2global.end(), i);
//                     int localIdx = it-mapFEM2global.begin();

//                     eparttmp[i] = epartMetis[localIdx];

//                 }
//                 else
//                 {
//                     it=find(mapMM2global.begin(), mapMM2global.end(), i);
//                     int localIdx = it-mapMM2global.begin();

//                     eparttmp[i] = maxProcNum+1+epartMetisMM[localIdx];
//                 }
//             }

//             for(int i =0; i<nNod;i++)
//             {
//                 bool exists = std::find(std::begin(nodesFEM), std::end(nodesFEM), i) != std::end(nodesFEM);

//                 if(exists)
//                 {
//                     it=find(nodesFEM.begin(), nodesFEM.end(), i);
//                     int localIdx = it-nodesFEM.begin();
//                     nparttmp[i] = npartMetis[localIdx];

//                 }
//                 else
//                 {
//                     it=find(nodesMM.begin(), nodesMM.end(), i);
//                     int localIdx = it-nodesMM.begin();
//                     nparttmp[i] = maxProcNum+1+npartMetisMM[localIdx];
//                 }
//             }
//         }
//         else
//         {
//             METIS_PartMeshNodal(&nElemMetis,&nNodMetis,mapmetisMetis,elemetisMetis,NULL,NULL,&nprcsMetis,NULL,NULL,&nMetis,epartMetis,npartMetis);
//             for(int i=0; i<nNod; i++)
//             {
//                 nparttmp[i] = npartMetis[i];
//             }

//             for(int i=0; i<nElem; i++)
//             {
//                 eparttmp[i] = epartMetis[i];

//             }

//         }
//     }

//     MPI_Barrier(MPI_COMM_WORLD); //Wait all the processors to reach this point
//     MPI_Bcast(&maxProcNum,1,MPI_INT,0,MPI_COMM_WORLD);

//     MPI_Bcast(nparttmp,nNod,MPI_INT,0,MPI_COMM_WORLD);
//     MPI_Bcast(eparttmp,nElem,MPI_INT,0,MPI_COMM_WORLD);

//     if(rank > maxProcNum)
//     {
//         flagProc4MM = true;
//     }

//     vector<int> npart2(nparttmp, nparttmp + nNod);
//     vector<int> epart2(eparttmp, eparttmp + nElem);

//     npart = npart2;
//     epart = epart2;

//     npartShare.resize(nNod, 0);
//     //Building the elements and nodes mapping vector and  the processor ownership vectors.
//     getownership(argc,argv,npart,epart,nodLocal,elemLocal,npartShare);

// }


//Dongli's changes with the option for the pre defined mesh partitions in fem or mm
extern void parallelManagement(int &nNodGlobal, int &nNod, int &nElem, vector<int> &epart,
                               vector<int> &npartShare, vector<int> &nodLocal, vector<int> &elemLocal,
                               vector<int> &npart, bool &flagProc4MM) {
    int nprcs = 1; // MPI number of processor and total number of processors, respectively
    int nodesPerElement = 0;
    nprcs = 0;
    int rank = 0;
    //Recover rank and number of processes information
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);//CHKERRQ(ierr);
    MPI_Comm_size(MPI_COMM_WORLD, &nprcs);//CHKERRQ(ierr);
    int nNodeFEM = 0;
    int nElemFEM = 0;
    int nNodeMM = 0;
    vector<int> elementsFEM;
    vector<int> nodesFEM;
    vector<int> nodesMM;
    bool flagCombined = false;
    int procNumMM = 1;
    bool flagMETIS = true;

    if (rank == 0) {
        INFO("Running simulation on %i cores", nprcs);
        INFO("Partitioning mesh ...");
        modelInfo(nNod, nElem, nodesPerElement, elementsFEM, nodesFEM, nodesMM, flagMETIS);
        nNodGlobal = nNod;
        nNodeFEM = nodesFEM.size();
        nElemFEM = elementsFEM.size();
        nNodeMM = nodesMM.size();
        /// This assigment of partition numbers for MM should be used when the partitioning application is able to generate convex partitions
        //procNumMM = int(nprcs * nElemMM * GPsPerMM / (nElemMM * GPsPerMM + nElemFEM * GPsPerFEM));//XXX
        /// Since METIS may generate non-convex partitions, we manually set the number of partitions to one for MM (this single partition should be convex)
        procNumMM = 1;
    }

    //Casting the information from the root process
    MPI_Barrier(MPI_COMM_WORLD); //Wait all the processors to reach this point
    MPI_Bcast(&nNod, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nNodGlobal, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nElem, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nodesPerElement, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nNodeFEM, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nElemFEM, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nNodeMM, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&procNumMM, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&flagMETIS, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (nNodeFEM == 0) {
        flagProc4MM = true;
    }

    if (flagMETIS) {
#ifdef METIS
        //Creating the METIS variables
        idx_t nElemMetis = nElemFEM;//nElem;                               // Number of elements in the global model
        idx_t nNodMetis = nNodeFEM; //nNod;                                 // Number of nodes in the global model
        idx_t nElemMetisMM = nElem - nElemFEM;
        idx_t nNodMetisMM = nNodeMM;
        idx_t nprcsMetis = nprcs - procNumMM;                              // default for combined FEM MM
        idx_t nprcsMetisMM = procNumMM;


        vector<int> mapMM2global(nElem - nElemFEM, 0);

        int maxProcNum = 0;//nprcsMetis-procNumMM;

        if (nNodeFEM == nNod) {
            nprcsMetis = nprcs;                              // FEM only
        } else if (nNodeFEM == 0) //MM only
        {
            // for MM

            nElemMetis = nElem - nElemFEM;
            nNodMetis = nNodeMM;
            nprcsMetis = nprcs;                              // FEM only
//            flagProc4MM = true;
        } else {
            flagCombined = true;
            nElemMetisMM = nElem - nElemFEM;
            nNodMetisMM = nNodeMM;
        }

        // Partitioning of the mesh by Metis
        vector<idx_t> npartMetis(nNodMetis);                         // Processor mapping of the nodes
        //assert(npartMetis != NULL);

        idx_t nMetis = 0;                                        // Number of Parts in Metis
        vector<idx_t> epartMetis(nElemMetis);                       // Processor mappingof th elements
        //assert(epartMetis != NULL);
        vector<int> mapFEM2global(nElemFEM, 0);
        vector<idx_t> nparttmp(nNod);
        vector<idx_t> eparttmp(nElem);

        // for MM
        MPI_Bcast(&nNodeMM, 1, MPI_INT, 0, MPI_COMM_WORLD);

        vector<idx_t> npartMetisMM(nNodMetisMM);                         // Processor mapping of the nodes
        //assert(npartMetisMM != NULL);

        idx_t nMetisMM = 0;                                        // Number of Parts in Metis
        vector<idx_t> epartMetisMM(nElemMetisMM);                       // Processor mappingof th elements
        //assert(epartMetisMM != NULL);

        vector<int>::iterator it;

        if (rank == 0) {
            vector<idx_t> elemetisMetis(nodesPerElement * nElemMetis); // Vector of nodes of all the elements
            vector<idx_t> mapmetisMetis(nElemMetis + 1); // Pre-metis vector distribution of elements in elemetisMetis
            vector<idx_t> elemetisMetisMM(nodesPerElement * nElemMetisMM);  // Vector of nodes of all the elements
            vector<idx_t> mapmetisMetisMM(nElemMetisMM + 1);  // Pre-metis vector distribution of elements in elemetisMetis
            
            if (flagCombined)
            {
                metisInput(&elemetisMetis[0], &mapmetisMetis[0], flagCombined, elementsFEM, mapFEM2global, nNodeFEM,
                       mapMM2global, &mapmetisMetisMM[0], &elemetisMetisMM[0], nNodeMM);
            }
            else
            {
                metisInput(&elemetisMetis[0], &mapmetisMetis[0], flagCombined, elementsFEM, mapFEM2global, nNodeFEM,
                       mapMM2global, NULL, NULL, nNodeMM);
            }
            
            if (flagCombined) {
                if (nprcsMetis == 1){
                    INFO("Number of partitions for FEM nodes and elements is set to 1 (METIS is not used)");
                    for (int i=0; i< nNodMetis; i++)
                        npartMetis[i] = 0;
                    for (int i=0; i< nElemMetis; i++)
                        epartMetis[i] = 0;
                    }
               else
                    METIS_PartMeshNodal(&nElemMetis, &nNodMetis, &mapmetisMetis[0], &elemetisMetis[0], NULL, NULL, &nprcsMetis,
                                        NULL, NULL, &nMetis, &epartMetis[0], &npartMetis[0]);
               
               if (nprcsMetisMM == 1){
                   INFO("Number of partitions for MM nodes and elements is set to 1 (METIS is not used)");
                   for (int i=0; i< nNodMetisMM; i++)
                        npartMetisMM[i] = 0;
                   for (int i=0; i< nElemMetisMM; i++)
                        epartMetisMM[i] = 0;
                   }
               else
                   METIS_PartMeshNodal(&nElemMetisMM, &nNodMetisMM, &mapmetisMetisMM[0], &elemetisMetisMM[0], NULL, NULL,
                                    &nprcsMetisMM, NULL, NULL, &nMetisMM, &epartMetisMM[0], &npartMetisMM[0]);
                

                vector<int> tmpConv(&epartMetis[0], &epartMetis[0] + nElemMetis);
                for (int k = 0; k < tmpConv.size(); k++) {
                    if (tmpConv[k] > maxProcNum) {
                        maxProcNum = tmpConv[k];
                    }
                }

                for (int i = 0; i < nElem; i++) {
                    bool exists =
                            std::find(std::begin(mapFEM2global), std::end(mapFEM2global), i) != std::end(mapFEM2global);

                    if (exists) {
                        it = find(mapFEM2global.begin(), mapFEM2global.end(), i);
                        int localIdx = it - mapFEM2global.begin();

                        eparttmp[i] = epartMetis[localIdx];

                    } else {
                        it = find(mapMM2global.begin(), mapMM2global.end(), i);
                        int localIdx = it - mapMM2global.begin();

                        eparttmp[i] = maxProcNum + 1 + epartMetisMM[localIdx];
                    }
                }

                for (int i = 0; i < nNod; i++) {
                    bool exists = std::find(std::begin(nodesFEM), std::end(nodesFEM), i) != std::end(nodesFEM);

                    if (exists) {
                        it = find(nodesFEM.begin(), nodesFEM.end(), i);
                        int localIdx = it - nodesFEM.begin();
                        nparttmp[i] = npartMetis[localIdx];

                    } else {
                        it = find(nodesMM.begin(), nodesMM.end(), i);
                        int localIdx = it - nodesMM.begin();
                        nparttmp[i] = maxProcNum + 1 + npartMetisMM[localIdx];
                    }
                }
            } else {

                METIS_PartMeshNodal(&nElemMetis, &nNodMetis, &mapmetisMetis[0], &elemetisMetis[0], NULL, NULL, &nprcsMetis,
                                    NULL, NULL, &nMetis, &epartMetis[0], &npartMetis[0]);
                for (int i = 0; i < nNod; i++) {
                    nparttmp[i] = npartMetis[i];
                }

                for (int i = 0; i < nElem; i++) {
                    eparttmp[i] = epartMetis[i];

                }

            }
        }

        MPI_Barrier(MPI_COMM_WORLD); //Wait all the processors to reach this point
        MPI_Bcast(&maxProcNum, 1, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Bcast(&nparttmp[0], nNod, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&eparttmp[0], nElem, MPI_INT, 0, MPI_COMM_WORLD);

        if (rank > maxProcNum) {
            flagProc4MM = true;
        }

        vector<int> npart2(&nparttmp[0], &nparttmp[0] + nNod);
        vector<int> epart2(&eparttmp[0], &eparttmp[0] + nElem);

        npart = npart2;
        epart = epart2;
#else
        ERROR("The code must be compiled with METIS to perform mesh partition at run-time");
        exit(-1);
#endif 
    } 
    else {
        std::cout << "predefined mesh" << endl;
        std::cout << "nprcs " << nprcs << endl;
        std::cout << "nNod " << nNod << endl;
        std::cout << "nElem " << nElem << endl;

        vector<int> flagFMMPI(nprcs);
        vector<int> npartMPI(nNod);
        vector<int> epartMPI(nElem);

        vector<int> arr_flagFM0(nprcs, 0);
        vector<int> npart0(nNod, 0);
        vector<int> epart0(nElem, 0);

        if (rank == 0) {
            partitionAssignment(npart0, epart0, arr_flagFM0, elementsFEM);
            std::cout << "partition assignment" << endl;
            std::cout << "npartMPI" << endl;
            for (int i = 0; i < nNod; i++) {
                npartMPI[i] = npart0[i];
            }
            std::cout << "epartMNPI" << endl;
            for (int i = 0; i < nElem; i++) {
                epartMPI[i] = epart0[i];
            }
            std::cout << "FM flag" << endl;
            for (int i = 0; i < nprcs; i++) {
                flagFMMPI[i] = arr_flagFM0[i];
            }
        }
        std::cout << "finish rank 0 read" << endl;
        MPI_Barrier(MPI_COMM_WORLD); //Wait all the processors to reach this point
        MPI_Bcast(&npartMPI[0], nNod, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&epartMPI[0], nElem, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&flagFMMPI[0], nprcs, MPI_INT, 0, MPI_COMM_WORLD);
        std::cout << "MPI cast " << rank << endl;

        if (flagFMMPI[rank] == 2) {
            flagProc4MM = true;
        }

        vector<int> npart2(&npartMPI[0], &npartMPI[0] + nNod);
        vector<int> epart2(&epartMPI[0], &epartMPI[0] + nElem);

        npart = npart2;
        epart = epart2;
    }
    std::cout << "finishing partition assign" << rank << endl;
    npartShare.resize(nNod, 0);
    for (int ii =0; ii < nNod; ii++)
    {
        if (npart[ii]==rank)
        {
            npartShare[ii] = 1;
        }
    }
    //Building the elements and nodes mapping vector and  the processor ownership vectors.
    getownership(npart, epart, nodLocal, elemLocal, npartShare);
    std::cout << "ownership rank " << rank << endl;

}

/*! \brief It assigns elements to processors according to input file //Sylvin, to correct in/out
  @param[in] argc Command line size
  @param[in] argv Command line data
  @param[out] nNod Number of global nodes
  @param[out] nNod Number of nodes
  @param[out] nElem number of Elements
  @param[out] epart Processor ownership mapping of the elements in the whole model epart[element i]=rank (empty for sequential simulations)
  @param[out] npartShare Processor ownership of the nodes in the whole model npartShare[node i]=1 if the nodes is mine or shared with other processor, 0 otherwise (empty for sequential simulations)
  @param[out] nodLocal Mapping of local nodes to global nodes numbering nodlocal[local node index]= global node index (empty for sequential simulations)
  @param[out] elemLocal Mapping of local elements to global elements numbering -> elemlocal[local element index] = global element index (empty for sequential simulations)
  @param[out] npart Processor ownership mapping of the nodes in the whole model npart[node i]=rank (empty for sequential simulations)
*/

void partitionAssignment(vector<int> &npart, vector<int> &epart, vector<int> &arr_flagFM,
                         vector<int> elemFEM) {
    ifstream ifs;
    string str = GeneralOptions::getInputFilePath();

    ifs.open(str.c_str(), ifstream::in);
    if (ifs.is_open()) {
        readPredef_part(ifs, npart, epart);
    } else {
        ERROR("The input file cannot be opened");
        exit(-1);
    }
   
    for (int i = 0; i < arr_flagFM.size(); i++) {
        int eleNum = std::find(epart.begin(), epart.end(), i) - epart.begin() + 1;
        bool exists = std::find(std::begin(elemFEM), std::end(elemFEM), eleNum) != std::end(elemFEM);
        if (exists) {
            arr_flagFM[i] = 1;//FEM
        } else {
            arr_flagFM[i] = 2;//MM
        }

    }
}

#endif


