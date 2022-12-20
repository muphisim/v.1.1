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

#ifndef _dataExtraction_H_
#define _dataExtraction_H_

#include <vector>
#include <string>
#include <map>
#include "classGPs.h"


struct ExtractData
{
    enum NodeDataType {NoneNodeType=0, Unknown, ExternalForce, InternalForce};
    enum ElementDataType {NoneEleType=0, FXX, FXY, FXZ, FYX, FYY, FYZ, FZX, FZY, FZZ,
                                  PXX, PXY, PXZ, PYX, PYY, PYZ, PZX, PZY, PZZ, SVM, KVM,
                                  INVAR1,INVAR2, INVAR3,INVAR4,INVAR5,INVAR6,
                                  DEFO_ENERGY, EXTERNAL_ENERGY, KIN_ENERGY, VOLUME
                                  };
    enum Operation {Undefined=0, Mean,Min,Max,Sum,Rough};
    static std::string to_string(NodeDataType data);
    static NodeDataType from_string_nodedatatype(std::string name);
    static std::string to_string(ElementDataType data);
    static ElementDataType from_string_elementdatatype(std::string name);
    static std::string to_string(Operation opt);
    static Operation from_string_operation(std::string opt);
};

class classDataBase
{
protected:
    std::string _name;
    FILE* _pFile;
    ExtractData::Operation  _operation;
    double _lastTime; 
        
public:
    classDataBase(std::string prefix, std::string opt);
    virtual ~classDataBase();
    
    double getLastTime() {return _lastTime;}
    std::string getName() const {return _name;};
    ExtractData::Operation getOperation() const {return _operation;};
    virtual void openFile();
    void closeFile();
};


class classNodeData : public classDataBase
{
protected:
    std::vector<int> _nodes;    
    ExtractData::NodeDataType _dataType;
    int _comp;
    
public:
    classNodeData(std::string prefix, const std::vector<int>& n, std::string type, std::string opt, int comp);
    virtual ~classNodeData();
    int getComp() const {return _comp;}
    ExtractData::NodeDataType getDataType() const {return _dataType;};
    const std::vector<int>& getNodes() const {return _nodes;}; 
    double getValue(const std::vector<double>& vals) const;
    void writeData(double time, const std::vector<double>& vals);
    
    virtual void openFile();
};

class classElementData : public classDataBase
{
protected:
    std::vector<int> _elements;    
    ExtractData::ElementDataType _dataType;
    bool _isInitialised;
    std::vector<classGPs*> _allGPs;
    std::vector<double> _weights;
    double _volume;
    bool _isAll;
    
public:
    classElementData(std::string type,  std::string opt);
    classElementData(std::string name, const std::vector<int>& eles, std::string type,  std::string opt);
    
    virtual ~classElementData();
    bool isInitialised() const {return _isInitialised;}
    std::vector<classGPs*>& getAllGPs() {return _allGPs;};
    ExtractData::ElementDataType getDataType() const {return _dataType;};
    const std::vector<int>& getElements() const {return _elements;}; 
    double getValue(const std::vector<double>& vals) const;
    void writeData(double time, const std::vector<double>& vals);
    void writeData(double time, double val);
    void initialise(std::vector<classGPs *> &GPs);
};


class classExtractData
{
protected:
    std::map<std::string, classNodeData*> _nodeDataToSave;
    std::map<std::string, classElementData*> _elementDataToSave;

public:
    classExtractData();
    ~classExtractData();
    void addClassNodeDataToSave(classNodeData* data);
    void addClassElementDataToSave(classElementData* data);
    void clear();
    void writeDataFromSolution(std::vector<classGPs *> &GPs, 
                              std::vector<double>& uK, std::vector<double>& vK, std::vector<double>& vvK,  
                              std::vector<double> & intForce, std::vector<double> & extForce,
                              int numNodes, double time);
    
public:
    static bool available;
    static classExtractData instance;
    static double getValueAtGP(classGPs* gpt, ExtractData::ElementDataType dataType);
    static double defoEnergy(std::vector<classGPs *> &GPs, ExtractData::Operation opt);
    static double kineticEnergy(std::vector<classGPs *> &GPs, std::vector<double>& vK, int numNodes, ExtractData::Operation opt);
    static double volume(std::vector<classGPs *> &GPs,  ExtractData::Operation opt);
    static void addData(classNodeData* data);
    static void addData(classElementData* data);
    static void extractField(std::vector<classGPs *> &GPs, 
                              std::vector<double>& uK, std::vector<double>& vK, std::vector<double>& vvK,  
                              std::vector<double> & intForce, std::vector<double> & extForce,
                              int numNodes, double time);
    static void finalise();
};

#endif //_dataExtraction_H_
