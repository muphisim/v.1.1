//
//
// File authors: see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file dataExtraction.h
  \brief This file contains all functions related the extract the data to the csv file
*/

#include "dataExtraction.h"
#include "commandLine.h"
#include "solvers.h"
#include <numeric>

/*! \brief This function is used to convert from NodeDataType to string
  @param[in] data NodeDataType type among Unknown, ExternalForce, and InternalForce
*/
std::string ExtractData::to_string(ExtractData::NodeDataType data)
{
    if (data == Unknown) return "Unknown";
    else if (data == ExternalForce) return "ExternalForce";
    else if (data == InternalForce) return "InternalForce";
    else return "NoneNodeType";
}
/*! \brief This function is used to convert a string to NodeDataType
  @param[in] name string name of a NodeDataType
*/
ExtractData::NodeDataType ExtractData::from_string_nodedatatype(std::string name)
{
    if (name == "Unknown") return Unknown;
    else if (name == "ExternalForce") return ExternalForce;
    else if (name == "InternalForce") return InternalForce;
    else
    {
        ERROR("node datatype %s has not defined",name.c_str());
        return NoneNodeType;
    }
}
/*! \brief This function is used to convert from ElementDataType to string
  @param[in] data ElementDataType 
*/
std::string ExtractData::to_string(ExtractData::ElementDataType data)
{
    if (data == ExtractData::FXX) return "FXX";
    else if (data == ExtractData::FXY) return "FXY";
    else if (data == ExtractData::FXZ) return "FXZ";
    else if (data == ExtractData::FYX) return "FYX";
    else if (data == ExtractData::FYY) return "FYY";
    else if (data == ExtractData::FYZ) return "FYZ";
    else if (data == ExtractData::FZX) return "FZX";
    else if (data == ExtractData::FZY) return "FZY";
    else if (data == ExtractData::FZZ) return "FZZ";
    else if (data == ExtractData::PXX) return "PXX";
    else if (data == ExtractData::PXY) return "PXY";
    else if (data == ExtractData::PXZ) return "PXZ";
    else if (data == ExtractData::PYX) return "PYX";
    else if (data == ExtractData::PYY) return "PYY";
    else if (data == ExtractData::PYZ) return "PYZ";
    else if (data == ExtractData::PZX) return "PZX";
    else if (data == ExtractData::PZY) return "PZY";
    else if (data == ExtractData::PZZ) return "PZZ";
    else if (data == ExtractData::SVM) return "SVM";
    else if (data == ExtractData::KVM) return "KVM";
    else if (data == ExtractData::INVAR1) return "INVAR1";
    else if (data == ExtractData::INVAR2) return "INVAR2";
    else if (data == ExtractData::INVAR3) return "INVAR3";
    else if (data == ExtractData::INVAR4) return "INVAR4";
    else if (data == ExtractData::INVAR5) return "INVAR5";
    else if (data == ExtractData::INVAR6) return "INVAR6";
    else if (data == ExtractData::DEFO_ENERGY) return "DEFO_ENERGY";
    else if (data == ExtractData::EXTERNAL_ENERGY) return "EXTERNAL_ENERGY";
    else if (data == ExtractData::KIN_ENERGY) return "KIN_ENERGY";
    else if (data == ExtractData::VOLUME) return "VOLUME";
    else return "NoneEleType";
}
/*! \brief This function is used to convert from a string to ElementDataType
  @param[in] name string name of a ElementDataType
*/
ExtractData::ElementDataType ExtractData::from_string_elementdatatype(std::string name)
{
    if (name == "FXX") return ExtractData::FXX;
    else if (name == "FXY") return ExtractData::FXY;
    else if (name == "FXZ") return ExtractData::FXZ;
    else if (name == "FYX") return ExtractData::FYX;
    else if (name == "FYY") return ExtractData::FYY;
    else if (name == "FYZ") return ExtractData::FYZ;
    else if (name == "FZX") return ExtractData::FZX;
    else if (name == "FZY") return ExtractData::FZY;
    else if (name == "FZZ") return ExtractData::FZZ;
    else if (name == "PXX") return ExtractData::PXX;
    else if (name == "PXY") return ExtractData::PXY;
    else if (name == "PXZ") return ExtractData::PXZ;
    else if (name == "PYX") return ExtractData::PYX;
    else if (name == "PYY") return ExtractData::PYY;
    else if (name == "PYZ") return ExtractData::PYZ;
    else if (name == "PZX") return ExtractData::PZX;
    else if (name == "PZY") return ExtractData::PZY;
    else if (name == "PZZ") return ExtractData::PZZ;
    else if (name == "SVM") return ExtractData::SVM;
    else if (name == "KVM") return ExtractData::KVM;
    else if (name == "INVAR1") return ExtractData::INVAR1;
    else if (name == "INVAR2") return ExtractData::INVAR2;
    else if (name == "INVAR3") return ExtractData::INVAR3;
    else if (name == "INVAR4") return ExtractData::INVAR4;
    else if (name == "INVAR5") return ExtractData::INVAR5;
    else if (name == "INVAR6") return ExtractData::INVAR6;
    else if (name == "DEFO_ENERGY") return ExtractData::DEFO_ENERGY;
    else if (name == "EXTERNAL_ENERGY") return ExtractData::EXTERNAL_ENERGY;
    else if (name == "KIN_ENERGY") return ExtractData::KIN_ENERGY;
    else if (name == "VOLUME") return ExtractData::VOLUME;
    else
    {
        ERROR("element datatype %s has not defined",name.c_str());
        return NoneEleType;
    }
}
/*! \brief The function is used to convert from Operation to string
  @param[in] opt Operation 
*/
std::string ExtractData::to_string(ExtractData::Operation opt)
{
    if (opt == Mean) return "Mean";
    else if (opt == Min) return "Min";
    else if (opt == Max) return "Max";
    else if (opt == Sum) return "Sum";
    else if (opt == Rough) return "Rough";
    else return "Undefined";
};
/*! \brief The function is used to convert a string from Operation
  @param[in] name string name of a Operation
*/
ExtractData::Operation ExtractData::from_string_operation(std::string opt)
{
    if (opt == "Mean") return Mean;
    else if (opt == "Min") return Min;
    else if (opt == "Max") return Max;
    else if (opt == "Sum") return Sum;
    else if (opt == "Rough") return Rough;
    else
    {
        ERROR("operation %s has not defined",opt.c_str());
        return Undefined;
    }
};

/*! \brief Constructor of the class
  @param[in] name string name of classDataBase
  @param[in] name string name of a Operation
*/
classDataBase::classDataBase(std::string name, std::string opt):
    _name(name),_pFile(NULL),_operation(ExtractData::from_string_operation(opt)), _lastTime(0.)
{
}

classDataBase::~classDataBase()
{
    closeFile();
}
    
void classDataBase::openFile()
{
    if (_pFile !=NULL) closeFile();
    std::string fname = "";
    int nprcs = 1;
    int rank = 0;
#ifdef PARALLEL 
    MPI_Comm_size(PETSC_COMM_WORLD, &nprcs);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
#endif    
    if (nprcs==1)
    {
        fname = GeneralOptions::getOutputDirectoryPath()+"/"+ getName()+ ".csv";
    }
    else
    {
        fname = GeneralOptions::getOutputDirectoryPath()+"/"+getName()+"_rank"+std::to_string(rank)+".csv";
    }
    _pFile = fopen(fname.c_str(),"w");
    fprintf(_pFile,"Time,Value\n");
};

void classDataBase::closeFile()
{
    if (_pFile !=NULL)
    {
        fclose(_pFile);
        _pFile = NULL;
    };
};

/*! \brief Constructor of the class
  @param[in] prefix string used as a first part of the file name
  @param[in] n list of nodes
  @param[in] type string of a NodeDataType
  @param[in] opt  string of a Operation
  @param[in] comp  component to save
*/
classNodeData::classNodeData(std::string prefix, const std::vector<int>& n, std::string type, std::string opt, int comp):
    classDataBase(prefix,opt),_nodes(n), _dataType(ExtractData::from_string_nodedatatype(type)), _comp(comp)
{
    _name = prefix+"_"+ExtractData::to_string(_dataType)+std::to_string(_comp)+"_"+ExtractData::to_string(_operation);
}

classNodeData::~classNodeData()
{
    
};

void classNodeData::openFile()
{
    if (_operation == ExtractData::Rough)
    {
        if (_pFile !=NULL) closeFile();
        std::string fname = "";
        int nprcs = 1;
        int rank = 0;
    #ifdef PARALLEL 
        MPI_Comm_size(PETSC_COMM_WORLD, &nprcs);
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    #endif    
        if (nprcs==1)
        {
            fname = GeneralOptions::getOutputDirectoryPath()+"/"+ getName()+ ".csv";
        }
        else
        {
            fname = GeneralOptions::getOutputDirectoryPath()+"/"+getName()+"_rank"+std::to_string(rank)+".csv";
        }
        _pFile = fopen(fname.c_str(),"w");
        fprintf(_pFile,"Time");
        for (int i=0; i< _nodes.size(); i++)
        {
            string name = ",Node"+std::to_string(_nodes[i]+1)+"Comp"+std::to_string(_comp);
            fprintf(_pFile,name.c_str());
        }
        fprintf(_pFile,"\n");
    }
    else
    {
        classDataBase::openFile();
    }
};

/*! \brief This function is used to estimate the output values from a list of values obtained with the nodes 
  @param[in] vals value of all nodes
*/
double classNodeData::getValue(const std::vector<double>& vals) const
{
    double val = 0;
    int N = vals.size();
    if (_operation == ExtractData::Mean)
    {
        val = std::accumulate(vals.begin(), vals.end(), 0.)/(double)N;
    }
    else if (_operation == ExtractData::Max)
    {
        val = *std::max_element(vals.begin(), vals.end());
    }
    else if (_operation == ExtractData::Min)
    {
        val = *std::min_element(vals.begin(), vals.end());
    }
    else if (_operation == ExtractData::Sum)
    {
        val = std::accumulate(vals.begin(), vals.end(), 0.);
    }
    else
    {
        ERROR("operation %d has not been defined",_operation);
    }
    return val;
}



/*! \brief Write data to file
  @param[in] time current time
  @param[in] vals values for writing that must be used with classNodeData::getValue
*/
void classNodeData::writeData(double time, const std::vector<double>& vals)
{
    if (time <= _lastTime) return;
    _lastTime = time;
    if (_pFile ==NULL) openFile();
    if (_operation == ExtractData::Rough)
    {
        fprintf(_pFile,"%.16g",time);
        for (int i=0; i< _nodes.size(); i++)
        {
            fprintf(_pFile,",%.16g",vals[i]);
        }
        fprintf(_pFile,"\n");
        fflush(_pFile);
    } 
    else
    {
        double val = getValue(vals);
        fprintf(_pFile,"%.16g,%.16g\n",time,val);
        fflush(_pFile);        
    }
};

/*! \brief Constructor of the class, this used when the operation is performed over all avaiable elements
  @param[in] type string of a ElementDataType
  @param[in] opt  string of a Operation
*/
classElementData::classElementData(std::string type, std::string opt):
    classDataBase("All",opt), _elements(), _dataType(ExtractData::from_string_elementdatatype(type)), 
    _weights(), _volume(0), _isAll(true), _isInitialised(false)
{
    _name = "All_"+ExtractData::to_string(_dataType)+"_"+ExtractData::to_string(_operation);
    if (_operation == ExtractData::Rough)
    {
        ERROR("Rough data is not implemented for classElementData");
        exit(-1);
    }
}

/*! \brief Constructor of the class
  @param[in] prefix string used as a first part of the file name
  @param[in] neles list of elements
  @param[in] type string of a ElementDataType
  @param[in] opt  string of a Operation
*/
classElementData::classElementData(std::string prefix, const std::vector<int>& eles, std::string type, std::string opt):
    classDataBase(prefix,opt), _elements(eles), _dataType(ExtractData::from_string_elementdatatype(type)), 
    _weights(), _volume(0), _isAll(false), _isInitialised(false)
{
    _name = prefix+"_"+ExtractData::to_string(_dataType)+"_"+ExtractData::to_string(_operation);
}

classElementData::~classElementData()
{
   
};
/*! \brief compute the output values form the values over each GPs the the current list of GPs
  @param[in] vals value of all nodes
*/
double classElementData::getValue(const std::vector<double>& vals) const
{
    double val = 0;
    if (_operation == ExtractData::Mean)
    {
        std::vector<double> weightedVals(vals);
        for (int igp=0; igp < _allGPs.size(); igp++)
        {
            weightedVals[igp] *= _weights[igp];
        }
        val = std::accumulate(weightedVals.begin(), weightedVals.end(), 0.);
    }
    else if (_operation == ExtractData::Max)
    {
        val = *std::max_element(vals.begin(), vals.end());
    }
    else if (_operation == ExtractData::Min)
    {
        val = *std::min_element(vals.begin(), vals.end());
    }
    else if (_operation == ExtractData::Sum)
    {
        std::vector<double> weightedVals(vals);
        for (int igp=0; igp < _allGPs.size(); igp++)
        {
            weightedVals[igp] *= _weights[igp];
        }
        val = _volume*std::accumulate(weightedVals.begin(), weightedVals.end(), 0.);
    }
    else
    {
        ERROR("operation %d has not been defined",_operation);
    }
    return val;
}

/*! \brief write to file
  @param[in] time current time
  @param[in] val value for writing
*/
void classElementData::writeData(double time, double val)
{
    if (!_isInitialised)
    {
        ERROR("classElementData must be initialised first");
        return;
    }
    if (time <= _lastTime) return;
    _lastTime = time;
    if (_pFile ==NULL) openFile();
    fprintf(_pFile,"%.16g,%.16g\n",time,val);
    fflush(_pFile);
}


/*! \brief write the file
  @param[in] time current time
  @param[in] vals values for writing 
*/
void classElementData::writeData(double time, const std::vector<double>& vals)
{
    if (!_isInitialised)
    {
        ERROR("classElementData must be initialised first");
        return;
    }
    if (time <= _lastTime) return;
    _lastTime = time;
    if (_pFile ==NULL) openFile();
    double val = getValue(vals);
    fprintf(_pFile,"%.16g,%.16g\n",time,val);
    fflush(_pFile);
};

/*! \brief initilise the GPs of the current object
  @param[in] GPs all GPs in the model
*/
void classElementData::initialise(std::vector<classGPs *> &GPs)
{
    _isInitialised = true;
    _allGPs.clear();
    _volume = 0;
    _weights.clear();
    for (int i=0; i< GPs.size(); i++)
    {
        classGPs* gpt = GPs[i];
        int ide = gpt->getCell();
        if (_isAll || std::find(_elements.begin(),_elements.end(),ide) != _elements.end())
        {
             double dV = gpt->getJ()*gpt->getW();
            _allGPs.push_back(gpt);
            _weights.push_back(dV);
            _volume += dV;
        }
    }
    for (int i=0; i< _weights.size(); i++)
    {
        _weights[i] /= _volume;
    }
};

/*! \brief Constructor of the class
*/
classExtractData::classExtractData(){};
classExtractData::~classExtractData()
{
    clear();
}

/*! \brief clear all existing data
*/
void classExtractData::clear()
{
    for (std::map<std::string, classNodeData*>::iterator it =  _nodeDataToSave.begin(); it != _nodeDataToSave.end(); it++)
    {
        delete it->second;
    }
    _nodeDataToSave.clear();
    for (std::map<std::string, classElementData*>::iterator it =  _elementDataToSave.begin(); it != _elementDataToSave.end(); it++)
    {
        delete it->second;
    }
    _elementDataToSave.clear();
}

/*! \brief add classNodeData object
  @param[in] data classNodeData object
*/
void classExtractData::addClassNodeDataToSave(classNodeData* data)
{
    if (data == NULL) return;
    if (_nodeDataToSave.find(data->getName()) == _nodeDataToSave.end())
    {
        _nodeDataToSave[data->getName()] = data;
    }
};

/*! \brief add classElementData object
  @param[in] data classElementData object
*/
void classExtractData::addClassElementDataToSave(classElementData* data)
{
    if (data == NULL) return;
    if (_elementDataToSave.find(data->getName()) == _elementDataToSave.end())
    {
        _elementDataToSave[data->getName()] = data;
    }
}

bool classExtractData::available = false;
classExtractData classExtractData::instance;

void classExtractData::addData(classNodeData* data)
{
    classExtractData::available = true;
    classExtractData::instance.addClassNodeDataToSave(data);
};
void classExtractData::addData(classElementData* data)
{
    classExtractData::available = true;
    classExtractData::instance.addClassElementDataToSave(data);
};

/*! \brief This function allows to obtain the value at each GP by ElementDataType
  @param[in] gpt Gauss point
  @param[in] dataType ElementDataType
*/
double classExtractData::getValueAtGP(classGPs* gpt, ExtractData::ElementDataType dataType)
{
    if (gpt->getActivate()==0) return 0.;
    const std::vector<double>& F = gpt->getCurrentState().getDeformationGradient();
    const std::vector<double>& P = gpt->getCurrentState().getFirstPiolaKirchhoffStress();
    auto getSVM = [](classGPs& GP)
    {
        int ndim = GP.getCurrentState().getDimension();
        std::vector<double> Cauchy(ndim * ndim, 0);
        GP.getCurrentState().getCauchy(Cauchy);
        double VMS = 0;
        if (ndim == 3) {
            VMS = sqrt(0.5 * ((Cauchy[0] - Cauchy[4]) * (Cauchy[0] - Cauchy[4]) +
                              (Cauchy[4] - Cauchy[8]) * (Cauchy[4] - Cauchy[8]) +
                              (Cauchy[8] - Cauchy[0]) * (Cauchy[8] - Cauchy[0]) +
                              6 * (Cauchy[1] * Cauchy[1] + Cauchy[2] * Cauchy[2] + Cauchy[5] * Cauchy[5])));
        } else {
            VMS = sqrt(
                    Cauchy[0] * Cauchy[0] + Cauchy[3] * Cauchy[3] - Cauchy[0] * Cauchy[3] + 3 * Cauchy[1] * Cauchy[1]);
        }
        return VMS;
    };
    int ndim = gpt->getDimension();
    const vector<double>& invars = gpt->getCurrentState().getInternalVariables();
    if (dataType == ExtractData::FXX) return F[0 * ndim + 0];
    else if (dataType == ExtractData::FXY) return F[0 * ndim + 1];
    else if (dataType == ExtractData::FXZ) return (ndim ==3)? F[0 * ndim + 2]: 0.;
    else if (dataType == ExtractData::FYX) return F[1 * ndim + 0];
    else if (dataType == ExtractData::FYY) return F[1 * ndim + 1];
    else if (dataType == ExtractData::FYZ) return (ndim ==3)? F[1 * ndim + 2]: 0.;
    else if (dataType == ExtractData::FZX) return (ndim ==3)? F[2 * ndim + 0]: 0.;
    else if (dataType == ExtractData::FZY) return (ndim ==3)? F[2 * ndim + 1]: 0.;
    else if (dataType == ExtractData::FZZ) return (ndim ==3)? F[2 * ndim + 2]: 0.;
    else if (dataType == ExtractData::PXX) return P[0 * ndim + 0];
    else if (dataType == ExtractData::PXY) return P[0 * ndim + 1];
    else if (dataType == ExtractData::PXZ) return (ndim ==3)? P[0 * ndim + 2]: 0.;
    else if (dataType == ExtractData::PYX) return P[1 * ndim + 0];
    else if (dataType == ExtractData::PYY) return P[1 * ndim + 1];
    else if (dataType == ExtractData::PYZ) return (ndim ==3)? P[1 * ndim + 2]: 0.;
    else if (dataType == ExtractData::PZX) return (ndim ==3)? P[2 * ndim + 0]: 0.;
    else if (dataType == ExtractData::PZY) return (ndim ==3)? P[2 * ndim + 1]: 0.;
    else if (dataType == ExtractData::PZZ) return (ndim ==3)? P[2 * ndim + 2]: 0.;
    else if (dataType == ExtractData::SVM) return getSVM(*gpt);
    else if (dataType == ExtractData::KVM) return getSVM(*gpt)*determinantTensor(F, ndim);
    else if (dataType == ExtractData::INVAR1) return (invars.size()>1)?invars[0]:0;
    else if (dataType == ExtractData::INVAR2) return (invars.size()>2)?invars[1]:0;
    else if (dataType == ExtractData::INVAR3) return (invars.size()>3)?invars[2]:0;
    else if (dataType == ExtractData::INVAR4) return (invars.size()>4)?invars[3]:0;
    else if (dataType == ExtractData::INVAR5) return (invars.size()>5)?invars[4]:0;
    else if (dataType == ExtractData::INVAR6) return (invars.size()>6)?invars[5]:0;
    else
    {
        ERROR("value %d at GP has not been defined",dataType);
        return 0;
    }
};

/*! \brief compute deformation energy
  @param[in] GPs a list of Gauss points
  @param[in] opt operation 
*/
double classExtractData::defoEnergy(std::vector<classGPs *> &GPs, ExtractData::Operation opt)
{
    double ener = 0.;
    double volume = 0.;
    bool initVal = false;
    for (vector<classGPs *>::iterator ite = GPs.begin(); ite != GPs.end(); ++ite) {

        classGPs *GaussPoint = *ite;
        if ((GaussPoint->getActivate()) == 0)
            continue;
        
        vector<constitutiveModels *> GaussPointCons = GaussPoint->getConstitutiveManager().getConsModel();
        const vector<double>& FCur = GaussPoint->getCurrentState().getDeformationGradient();
        const vector<double>& intVarsCur = GaussPoint->getCurrentState().getInternalVariables();
        double wJ = GaussPoint->getWeightJ();
        double val = GaussPointCons[0]->defoEnergy(GaussPoint);

        double dV = wJ;
        volume += dV;
        if (opt == ExtractData::Sum || opt == ExtractData::Mean)
        {
            ener += val*dV;
        } 
        else if (opt == ExtractData::Max)
        {
            if (!initVal)
            {
                initVal = true;
                ener = val;
            }
            else if (val > ener)
            {
                ener = val;
            }
        }
        else if (opt == ExtractData::Min)  
        {
            if (!initVal)
            {
                initVal = true;
                ener = val;
            }
            else if (val < ener)
            {
                ener = val;
            }
        }
    }
    if (opt == ExtractData::Mean)
    {
        return ener/volume;
    }
    else
        return ener;
};

/*! \brief compute kinetic energy
  @param[in] GPs a list of Gauss points
  @param[in] vK vector of velocity
  @param[in] numNodes number of nodes
  @param[in] opt operation 
*/
double classExtractData::kineticEnergy(std::vector<classGPs *> &GPs, std::vector<double>& vK, int numNodes, ExtractData::Operation opt)
{
    double ener = 0.;
    double volume = 0.;
    bool initVal = false;
    for (vector<classGPs *>::iterator ite = GPs.begin(); ite != GPs.end(); ++ite) {

        classGPs *GaussPoint = *ite;
        if ((GaussPoint->getActivate()) == 0)
            continue;
        int dim = GaussPoint->getDimension();
        vector<constitutiveModels *> GaussPointCons = GaussPoint->getConstitutiveManager().getConsModel();
        const vector<double>& FCur = GaussPoint->getCurrentState().getDeformationGradient();
        const vector<double>& intVarsCur = GaussPoint->getCurrentState().getInternalVariables();
        double w = GaussPoint->getW();
        double Jacobian = GaussPoint->getJ();
        const std::vector<int>& neighbours = GaussPoint->getNeighbours();
        const vector<double>& phi = GaussPoint->getPhi();
        double rho = GaussPointCons[0]->getRho();
        std::vector<double> velocity(3,0.);
        for (int i=0; i< neighbours.size(); i++)
        {
            for (int j=0; j<dim; j++)
            {
                velocity[j] += phi[i]*vK[neighbours[i]+j*numNodes];
            }
            
        }
        double val = 0;
        for (int j=0; j< dim; j++)
        {
            val += 0.5*rho*velocity[j]*velocity[j];
        }
        double dV = w * Jacobian;
        volume += dV;
        if (opt == ExtractData::Sum || opt == ExtractData::Mean)
        {
            ener += val*dV;
        } 
        else if (opt == ExtractData::Max)
        {
            if (!initVal)
            {
                initVal = true;
                ener = val;
            }
            else if (val > ener)
            {
                ener = val;
            }
        }
        else if (opt == ExtractData::Min)  
        {
            if (!initVal)
            {
                initVal = true;
                ener = val;
            }
            else if (val < ener)
            {
                ener = val;
            }
        }
    }
    if (opt == ExtractData::Mean)
    {
        return ener/volume;
    }
    else
        return ener;
}

/*! \brief compute volume
  @param[in] GPs a list of Gauss points
  @param[in] vK vector of velocity
  @param[in] numNodes number of nodes
  @param[in] opt operation
*/
double classExtractData::volume(std::vector<classGPs *> &GPs, ExtractData::Operation opt)
{

    double volume = 0.;
    bool initVal = false;
    double J_mech = 0.;
    for (vector<classGPs *>::iterator ite = GPs.begin(); ite != GPs.end(); ++ite) {

        classGPs *GaussPoint = *ite;
        if ((GaussPoint->getActivate()) == 0)
            continue;
        int dim = GaussPoint->getDimension();
        vector<constitutiveModels *> GaussPointCons = GaussPoint->getConstitutiveManager().getConsModel();
        const vector<double>& FCur = GaussPoint->getCurrentState().getDeformationGradient();
        const vector<double>& intVarsCur = GaussPoint->getCurrentState().getInternalVariables();
        double w = GaussPoint->getW();
        double Jacobian = GaussPoint->getJ();
        const std::vector<int>& neighbours = GaussPoint->getNeighbours();
        const vector<double>& phi = GaussPoint->getPhi();
        J_mech = determinantTensor(FCur,dim);
        double dV = w * Jacobian * J_mech;
        if (opt == ExtractData::Sum || opt == ExtractData::Mean)
        {
            volume += dV;
        }
        else if (opt == ExtractData::Max)
        {
            if (!initVal)
            {
                initVal = true;
                volume = dV;
            }
            else if (dV > volume)
            {
                volume = dV;
            }
        }
        else if (opt == ExtractData::Min)
        {
            if (!initVal)
            {
                initVal = true;
                volume = dV;
            }
            else if (dV < volume)
            {
                volume = dV;
            }
        }
    }
    if (opt == ExtractData::Mean)
    {

    }
    else
        return volume;
}


/*! \brief write data to file
  @param[in] GPs a list of Gauss points
  @param[in] uK displacement
  @param[in] vK velocity
  @param[in] vvK extraDof
  @param[in] intForce internal force 
  @param[in] extForce external force
  @param[in] numNodes number of nodes
  @param[in] time current time 
*/
void classExtractData::writeDataFromSolution(std::vector<classGPs *> &GPs, 
                              std::vector<double>& uK, std::vector<double>& vK, std::vector<double>& vvK,   
                              std::vector<double> & intForce, std::vector<double> & extForce,
                              int numNodes, double time)
{
    for (std::map<std::string, classNodeData*>::iterator it =  _nodeDataToSave.begin(); it != _nodeDataToSave.end(); it++)
    {
        classNodeData* dataToSave = it->second;
        if (time > dataToSave->getLastTime())
        {
            const std::vector<int>& nodes = dataToSave->getNodes();
            int NSize = nodes.size();
            std::vector<double> values(NSize,0.);
            ExtractData::NodeDataType dataType =  dataToSave->getDataType();
            int numMechaDof = uK.size()/numNodes;
            int numExtraDof = vvK.size()/numNodes;
            for (int j=0; j <NSize; j++)
            {
                if (dataType == ExtractData::Unknown)
                {
                    int comp = dataToSave->getComp();
                    if (comp <numMechaDof)
                    {
                        values[j] = uK[nodes[j]+comp*numNodes];
                    }
                    else if (comp < numMechaDof+numExtraDof)
                    {
                        values[j] = vvK[nodes[j]+(comp-numMechaDof)*numNodes];
                    }
                    else
                    {
                        ERROR("component %d does not exist",comp);
                    }
                }
                else if (dataType == ExtractData::ExternalForce)
                {
                    int comp = dataToSave->getComp();
                    if (comp <numMechaDof)
                    {
                        values[j] = extForce[nodes[j]+comp*numNodes];
                    }
                    else
                    {
                        ERROR("component %d does not exist to get external force",comp);
                    }
                }
                else if (dataType == ExtractData::InternalForce)
                {
                    int comp = dataToSave->getComp();
                    if (comp <numMechaDof)
                    {
                        values[j] = intForce[nodes[j]+comp*numNodes];
                    }
                    else
                    {
                        ERROR("component %d does not exist to get internal force",comp);
                    }
                }
                
            }
            dataToSave->writeData(time,values);
        }
    };
        
    for (std::map<std::string, classElementData*>::iterator it =  _elementDataToSave.begin(); it != _elementDataToSave.end(); it++)
    {
        classElementData* dataToSave = it->second;
        if (!dataToSave->isInitialised())
        {
            dataToSave->initialise(GPs);
        }
        if (time > dataToSave->getLastTime())
        {
            ExtractData::ElementDataType dataType =  dataToSave->getDataType();
            if (dataType == ExtractData::DEFO_ENERGY)
            {
                double energ = classExtractData::defoEnergy(dataToSave->getAllGPs(),dataToSave->getOperation());
                dataToSave->writeData(time,energ);
            }
            else if (dataType == ExtractData::KIN_ENERGY) 
            {
                double energ = classExtractData::kineticEnergy(dataToSave->getAllGPs(),vK,numNodes,dataToSave->getOperation());
                dataToSave->writeData(time,energ);
            }
            else if (dataType == ExtractData::VOLUME)
            {
                double volume = classExtractData::volume(dataToSave->getAllGPs(),dataToSave->getOperation());
                dataToSave->writeData(time,volume);
            }
            else if (dataType == ExtractData::EXTERNAL_ENERGY)
            {
                double energ = 0;
                for (int i=0; i< uK.size(); i++)
                {
                    energ += uK[i]*extForce[i];
                }
                dataToSave->writeData(time,energ);
            }
            else
            {
                
                std::vector<classGPs*>& allGPs = dataToSave->getAllGPs();
                int NSize = allGPs.size();
                std::vector<double> values(NSize,0.);
                for (int j=0; j <NSize; j++)
                {
                    values[j] = classExtractData::getValueAtGP(allGPs[j],dataType);
                }
                dataToSave->writeData(time,values);
            }
        };
    };
};

/*! \brief write data to file
  @param[in] GPs a list of Gauss points
  @param[in] uK displacement
  @param[in] vK velocity
  @param[in] vvK extraDof
  @param[in] intForce internal force 
  @param[in] extForce external force
  @param[in] numNodes number of nodes
  @param[in] time current time
*/
void classExtractData::extractField(std::vector<classGPs *> &GPs, 
                               std::vector<double>& uK, std::vector<double>& vK, std::vector<double>& vvK,  
                              std::vector<double> & intForce, std::vector<double> & extForce,
                              int numNodes, double time)
{
    if (classExtractData::available == false) return; // do nothing
    classExtractData::instance.writeDataFromSolution(GPs,uK,vK,vvK,intForce,extForce,numNodes,time);
};

/*! \brief clear all existing data in this manager
*/
void classExtractData::finalise()
{
    classExtractData::available = false;
    classExtractData::instance.clear();
};
