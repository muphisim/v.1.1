//
//
// File author(s): see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _classNeumannBCs_H
#define _classNeumannBCs_H

#include "configuration.h"
#include "classSpatialPoint.h"
#include "constitutiveModels.h"

class classNeumannBCs {

protected:
    string _type;
public:
    classNeumannBCs();

    classNeumannBCs(string type);

    inline string getType() { return _type; }

    virtual ~classNeumannBCs() {};
};


/*! \brief  This class defines each Neumann node in the domain*/
class classNodalForces : public classNeumannBCs {

private:
    vector<classNeumann *> _listOfNodes;

public:
    classNodalForces(string type, vector<classNeumann *> listOfNodes);
    inline vector<classNeumann *>& getListOfNodes() { return _listOfNodes; }; // Common
    virtual ~classNodalForces() {
        /* for(int i=0;i<this->listOfNodes.size();i++){ */
        /*   delete this->listOfNodes[i]; */
        /* } */
    };
};

/*! \brief  This class defines each Neumann node in the domain*/
class classTractionInst : public classNeumannBCs {

private:
    vector<classElements *> _listOfElements; // They will be ndim-1 dimension
    double _Tx, _Ty, _Tz; // traction following 3 component
    double _t0, _tf; //
public:
    classTractionInst(string type, vector<classElements *> &listOfElements, double Tx, double Ty, double Tz, double t0, double tf);
    inline vector<classElements *>& getListOfElements() { return _listOfElements; }; // Common
    inline void getTraction(double timeRun, double dt, double& Tx, double& Ty, double& Tz) 
    {
        Tx =0;
        Ty =0;
        Tz =0;
        if (timeRun >= _t0 && timeRun <= _tf)
        {
            Tx = _Tx;
            Ty = _Ty;
            Tz = _Tz;
        } 
    }; // Common
    virtual ~classTractionInst() {

    };
};

/*! \brief  This class defines each Neumann node in the domain*/
class classTractionRamp : public classNeumannBCs {

private:
    vector<classElements *> _listOfElements; // They will be ndim-1 dimension
    double _Tx, _Ty, _Tz; // traction following 3 component
    double _t0, _tf; //
public:
    classTractionRamp(string type, vector<classElements *> &listOfElements, double Tx, double Ty, double Tz, double t0, double tf);
    inline vector<classElements *>& getListOfElements() { return _listOfElements; }; // Common
    inline void getTraction(double timeRun, double dt, double& Tx, double& Ty, double& Tz) 
    {
        Tx =0;
        Ty =0;
        Tz =0;
        if (timeRun >= _t0 && timeRun <= _tf)
        {
            Tx = _Tx / (_tf - _t0) * (timeRun - _t0);
            Ty = _Ty / (_tf - _t0) * (timeRun - _t0);
            Tz = _Tz / (_tf - _t0) * (timeRun - _t0);
        } 
    }; // Common
    virtual ~classTractionRamp() {

    };
};



/*! \brief  This class defines each Neumann node in the domain*/
class classStochasticPressureRamp : public classNeumannBCs {

private:
    vector<classElements *> _listOfElements; // They will be ndim-1 dimension
    double _t0, _tf;
    vector<double> _P, _C;
    int _nstoch;
public:
    classStochasticPressureRamp(string type, vector<classElements *> &listOfElements, double P, double P_var, double t0,
                                double tf, int nstoch, vector<double> C);

    inline vector<classElements *> getListOfElements() { return _listOfElements; }; // Common
    inline vector<double> getPressure(double timeRun, double dt) {
        vector<double> P(_nstoch,0);
        if (timeRun >= _t0 && timeRun <= _tf) {
            for (int i=0; i<_nstoch; i++) {
                P[i] = _P[i] / (_tf - _t0) * (timeRun - _t0);
            }
        } else {
            for (int i=0; i<_nstoch; i++) {
                P[i] = 0.;
            }
        }
        return P;
    }; // Common
    virtual ~classStochasticPressureRamp() {

    };
};

/*! \brief  This class defines each Neumann node in the domain*/
class classPressureRamp : public classNeumannBCs {

private:
    vector<classElements *> _listOfElements; // They will be ndim-1 dimension
    double _P, _t0, _tf;
public:
    classPressureRamp(string type, vector<classElements *> &listOfElements, double P, double t0, double tf);

    inline vector<classElements *>& getListOfElements() { return _listOfElements; }; // Common
    inline double getPressure(double timeRun, double dt) {
        double P = 0;
        if (timeRun >= _t0 && timeRun <= _tf) {
            P = _P / (_tf - _t0) * (timeRun - _t0);
        } else {
            P = 0.;
        }
        return P;
    }; // Common
    virtual ~classPressureRamp() {

    };
};

/*! \brief  This class defines each Neumann node in the domain*/
class classPressureInst : public classNeumannBCs {

private:
    vector<classElements *> _listOfElements; // They will be ndim-1 dimension
    double _P, _t0, _tf;
public:
    classPressureInst(string type, vector<classElements *> &listOfElements, double P, double t0, double tf);

    inline vector<classElements *>& getListOfElements() { return _listOfElements; }; // Common
    inline double getPressure(double timeRun, double dt) {
        double P = 0;
        if (timeRun >= _t0 && timeRun <= _tf) {
            P = _P;
        } else {
            P = 0.;
        }
        return P;
    }; // Common
    virtual ~classPressureInst() {

    };
};

/*! \brief  This class defines each Neumann element in the domain*/
class classCurrentInst : public classNeumannBCs {

private:
    vector<classElements *> _listOfElements; // They will be ndim-1 dimension
    double _INeumann, _t0, _tf;
    int _fieldIndex;
public:
    classCurrentInst(string type, vector<classElements *> &listOfElements, double P, double t0, double tf, int field=0);
    int getFieldIndex() const {return _fieldIndex;};
    inline vector<classElements *>& getListOfElements() { return _listOfElements; }; // Common
    inline double getCurrent(double timeRun, double dt) {
        double INeumann = 0;
        if (timeRun >= _t0 && timeRun <= _tf) {
            INeumann = _INeumann;
        } else {
            INeumann = 0.;
        }
        return INeumann;
    }; // Common
    virtual ~classCurrentInst() {

    };
};

/*! \brief  This class defines each Neumann node in the domain*/
class classHertzianRamp : public classNeumannBCs {

private:
    vector<classElements *> _listOfElements; // They will be ndim-1 dimension
    double _P0, _Q0, _a, _mu, _c,  _t0, _tf;
    vector<double> _center;
public:
    classHertzianRamp(string type, vector<classElements *> &listOfElements, double P, double Q, double mu, double a, vector<double> center, double t0, double tf);

    inline vector<classElements *>& getListOfElements() { return _listOfElements; }; // Common
    inline double getNormalPressure(double timeRun, double dt) {
        double P = 0;
        if (timeRun >= _t0 && timeRun <= _tf) {
            P = _P0 / (_tf - _t0) * (timeRun - _t0);
        } else {
            P = 0.;
        }
        return P;
    };
    inline double getTangentialPressure(double timeRun, double dt) {
        double Q = 0;
        if (timeRun >= _t0 && timeRun <= _tf) {
            Q = _Q0 / (_tf - _t0) * (timeRun - _t0);
        } else {
            Q = 0.;
        }
        return Q;
    };
    inline double getRadiusStick() {
return _c;
    };
    inline double getRadius() {
        return _a;
    };
    inline double getFriction() {
        return _mu;
    };
    inline vector<double> getCenter() {
        return _center;
    };
    // Common
    ~classHertzianRamp() {

    };
};

/*! \brief  This class defines each Neumann node in the domain*/
class classStochasticHertzianRamp : public classNeumannBCs {

private:
    vector<classElements *> _listOfElements; // They will be ndim-1 dimension
    double _P0,  _a,  _t0, _tf;
    int _nstoch;
    vector<double> _center;
    vector<double> _mu;
    vector<double> _c;
    vector<double> _C;
    vector<double> _Q0;
public:

    classStochasticHertzianRamp(string type, vector<classElements *> &listOfElements, double P, double Q, double mu, double mu_var, double a, vector<double> center, double t0, double tf, int nstoch, vector<double> C, bool flagHaar, int resolution, string distribution);

    inline vector<classElements *>& getListOfElements() { return _listOfElements; }; // Common
    inline double getNormalPressure(double timeRun, double dt) {
        double P = 0;
        if (timeRun >= _t0 && timeRun <= _tf) {
            P = _P0 / (_tf - _t0) * (timeRun - _t0);
        } else {
            P = 0.;
        }
        return P;
    };
    inline vector<double> getTangentialPressure(double timeRun, double dt) {
        int size;
        size = _Q0.size();
        vector<double> Q(size,0);
        for(int i=0; i<size;i++){
            if (timeRun >= _t0 && timeRun <= _tf) {
                Q[i] = _Q0[i] / (_tf - _t0) * (timeRun - _t0);
            } else {
                ;
            }
        }
        return Q;
    };
    inline const vector<double>& getRadiusStick() const {
        return _c;
    };
    inline double getRadius() {
        return _a;
    };
    inline const vector<double>& getFriction() const{
        return _mu;
    };
    inline const vector<double>& getCenter() const {
        return _center;
    };
    // Common
    ~classStochasticHertzianRamp() {

    };
};

/*! \brief  This class defines the volumetric heat following exponential function
 * 
 * q(x,y,z, t) = qmax*exp(-RR)
 * xc = x0+ vx*(t-t0)
 * yc = y0+ vy*(t-t0)
 * zc = z0+ vz*(t-t0)
 * with RR = ((x-xc)/a)^2 ((y-yc)/b)^2+((z-zc)/c)^2
 * 
*/

class classVolHeatFluxGaussian : public classNeumannBCs {

private:
    double _qmax;
    double _x0, _y0, _z0; // laser initial position
    double _vx, _vy, _vz; // laser velocity
    double _a, _b, _c;
    double _t0, _tf;
    int _fieldIndex;
public:
    classVolHeatFluxGaussian(string type, double qmax, double x0, double y0, double z0,
                                            double vx, double vy, double vz,    
                                            double a, double b, double c,
                                            double t0, double tf, int field=0);
    int getFieldIndex() const {return _fieldIndex;};
    double get_t0() const {return _t0;};
    double get_tf() const {return _tf;};
    inline double getVolHeatFlux(double x, double y, double z, double timeRun, double dt)
    {
        if (timeRun >= _t0 && timeRun < _tf) 
        {
            double xc = _x0 + _vx*(timeRun-_t0);
            double yc = _y0 + _vy*(timeRun-_t0);
            double zc = _z0 + _vz*(timeRun-_t0);
            double RR = (x-xc)*(x-xc)/(_a*_a) + (y-yc)*(y-yc)/(_b*_b)+(z-zc)*(z-zc)/(_c*_c);
            return _qmax*exp(-RR);
        } else {
            return 0.;
        }
        
    }; // Common
    virtual ~classVolHeatFluxGaussian() {

    };
};




/*! \brief  This class defines each Neumann element in the domain*/
class classVolHeatFlux: public classNeumannBCs {

private:
    vector<classElements *> _listOfElements;
    double _r, _t0, _tf;
    int _fieldIndex;
public:
    classVolHeatFlux(string type, vector<classElements *> &listOfElements, double r, double t0, double tf, int field=0);
    int getFieldIndex() const {return _fieldIndex;};
    inline vector<classElements *>& getListOfElements() { return _listOfElements; };
    inline double getVolHeatFlux(double timeRun, double dt) {
        double r = 0;
        if (timeRun >= _t0 && timeRun < _tf) {
            r = _r;
        } else {
            r = 0.;
        }
        return r;
    }; // Common
    virtual ~classVolHeatFlux() {

    };
};

/*! \brief  This class defines each Heat equation Neumann element in the domain*/
class classSurfaceHeatFluxInst : public classNeumannBCs {

private:
    vector<classElements *> _listOfElements; // They will be ndim-1 dimension
    vector<double> _HeatNeumann;
    double _t0, _tf;
    int _fieldIndex;
public:

    classSurfaceHeatFluxInst(string type, vector<classElements *> &listOfElements, vector<double> HeatNeumann, double t0,
                      double tf, int field=0);

    inline vector<classElements *>& getListOfElements() { return _listOfElements; }; // Common
    int getFieldIndex() const {return _fieldIndex;};
    inline vector<double> getHeat(double timeRun, double dt) {
        vector<double> HeatNeumann;
        if (timeRun >= _t0 && timeRun <= _tf) {
            HeatNeumann = _HeatNeumann;
        } else {
            HeatNeumann =  { 0, 0, 0 };
        }
        return HeatNeumann;
    }; // Common
    virtual ~classSurfaceHeatFluxInst() {

    };
};

/*! \brief  This class defines Convection BC*/
class classConvection : public classNeumannBCs {

private:
    vector<classElements *> _listOfElements; // They will be ndim-1 dimension
    double _convectiveCoef, _sinkTemp, _t0, _tf;
    int _fieldIndex;
public:

    classConvection(string type, vector<classElements *> &listOfElements, double filmCoefficient, double sinkTemp,
                    double t0, double tf, int field = 0);

    inline vector<classElements *>& getListOfElements() { return _listOfElements; }; // Common
    int getFieldIndex() const {return _fieldIndex;};
    inline double getConvection(double timeRun, double dt, double theta) { //Surface convection
        double qConvection;
        if (timeRun >= _t0 && timeRun < _tf) {
            qConvection = _convectiveCoef * (theta - _sinkTemp);
        } else {
            qConvection = 0;
        }
        return qConvection;
    };

    inline double getConvCons(double timeRun, double dt) { //Function to get convection constant o film coefficient h
        double hConvection;;
        if (timeRun >= _t0 && timeRun < _tf) {
            hConvection = _convectiveCoef;
        } else {
            hConvection = 0;
        }
        return hConvection;
    };
    virtual ~classConvection() {

    };
};

/*! \brief  This class defines Radiation BC*/
class classRadiation : public classNeumannBCs {

private:
    vector<classElements *> _listOfElements; // They will be ndim-1 dimension
    double _radiationConstant, _sinkTemp, _zeroTemp, _t0, _tf;
    int _fieldIndex;
public:

    classRadiation(string type, vector<classElements *> &listOfElements, double radiationConstant, double sinkTemp,
                   double t0, double tf, int field =0);
    int getFieldIndex() const {return _fieldIndex;};
    inline vector<classElements *> getListOfElements() { return _listOfElements; }; // Common
    inline double getRadiation(double timeRun, double dt, double theta) {  //Surface radiation
        double qRadiation;;
        if (timeRun >= _t0 && timeRun < _tf) { //
            qRadiation = _radiationConstant * (pow(theta, 4.0) - pow(_sinkTemp, 4.0));
        } else {
            qRadiation = 0;
        }
        return qRadiation;
    };

    inline double getRadiaCons(double timeRun, double dt) { //Function to get radiation constant A
        double ARadiation;;
        if (timeRun >= _t0 && timeRun < _tf) {
            ARadiation = _radiationConstant;
        } else {
            ARadiation = 0;
        }
        return ARadiation;
    };
    virtual ~classRadiation() {

    };
};

class classGravity : public classNeumannBCs {
private:
    double _g;
public:
    inline double getGravity() {
        return _g;
    }
    classGravity(string type, double g);

    virtual ~classGravity() {

    };
};

/*! \brief  This class defines FSI BC*/
class classFSI : public classNeumannBCs {

private:
    vector<int> _boundaryNodes;
    vector<vector <double> > _boundaryCoordinates;
public:
    classFSI(string type, vector<int> &boundary_nodes, vector<vector <double> > boundary_coordinates);
    inline vector<int>& getListOfBoundaryNodes() { return _boundaryNodes;};
    inline vector<vector <double> >& getListOfBoundaryCoordinates() {return _boundaryCoordinates;}
    inline int getNumberOfBoundaryNodes() { return _boundaryNodes.size();};
    virtual ~classFSI() {

    };
};

#endif
