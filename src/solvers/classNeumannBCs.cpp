//
//
// File author(s): see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file classNeumannBCs.cpp
  \brief This file defines all functions needed to apply the Neumann boundary conditions
*/

#include "classNeumannBCs.h"
#include "maths.h"
/*! Constructor for all Neumann BCs
  @param[in] type Type of Neumann BCs
*/
classNeumannBCs::classNeumannBCs(string type) {
    _type = type;
};

/*! Constructor for classNodalForces
  @param[in] type Type of Neumann BCs
  @param[in] listOfNodes Nodes at which the forces are applied
*/
classNodalForces::classNodalForces(string type, vector<classNeumann *> listOfNodes) : classNeumannBCs(type) {
    _listOfNodes = listOfNodes;
};

/*! Constructor for Neumann Traction BCs Inst
  @param[in] type Type of Neumann BCs
  @param[in] listOfElements Nodes at which the forces are applied
  @param[in] P Pressure
  @param[in] t0 Init time to apply the pressure
  @param[in] tf Final time to apply the pressure
*/
classTractionInst::classTractionInst(string type, vector<classElements *> &listOfElements, 
                                     double Tx, double Ty, double Tz,
                                     double t0, double tf) : classNeumannBCs(type) {
    _Tx = Tx;
    _Ty = Ty;
    _Tz = Tz;
    _listOfElements = listOfElements;
    _tf = tf;
    _t0 = t0;

};


/*! Constructor for Neumann Traction BCs Ramp
  @param[in] type Type of Neumann BCs
  @param[in] listOfElements Nodes at which the forces are applied
  @param[in] P Pressure
  @param[in] t0 Init time to apply the pressure
  @param[in] tf Final time to apply the pressure
*/
classTractionRamp::classTractionRamp(string type, vector<classElements *> &listOfElements, 
                                     double Tx, double Ty, double Tz,
                                     double t0, double tf) : classNeumannBCs(type) {
    _Tx = Tx;
    _Ty = Ty;
    _Tz = Tz;
    _listOfElements = listOfElements;
    _tf = tf;
    _t0 = t0;

};


/*! Constructor for Neumann Pressure BCs Ramp
  @param[in] type Type of Neumann BCs
  @param[in] listOfElements Nodes at which the forces are applied
  @param[in] P Pressure
  @param[in] t0 Init time to apply the pressure
  @param[in] tf Final time to apply the pressure
*/
classStochasticPressureRamp::classStochasticPressureRamp(string type, vector<classElements *> &listOfElements, double P, double P_var, double t0,
                                     double tf, int nstoch, vector<double> C) : classNeumannBCs(type) {
    _nstoch = nstoch;
    _C = C;
    _listOfElements = listOfElements;
    _tf = tf;
    _t0 = t0;
    _P.resize(nstoch);
    _P[0] = P;
    _P[1] = P_var;
};

/*! Constructor for Neumann Pressure BCs Ramp
  @param[in] type Type of Neumann BCs
  @param[in] listOfElements Nodes at which the forces are applied
  @param[in] P Pressure
  @param[in] t0 Init time to apply the pressure
  @param[in] tf Final time to apply the pressure
*/
classPressureRamp::classPressureRamp(string type, vector<classElements *> &listOfElements, double P, double t0,
                                     double tf) : classNeumannBCs(type) {
    _P = P;
    _listOfElements = listOfElements;
    _tf = tf;
    _t0 = t0;

};

/*! Constructor for Hertzian Pressure BCs Ramp
  @param[in] type Type of Neumann BCs
  @param[in] listOfElements Nodes at which the forces are applied
  @param[in] P TOTAL normal force
  @param[in] Q TOTAL tangential force
  @param[in] mu Friction coefficient
  @param[in] a Radius of contact
  @param[in] center Coordinates of the center
  @param[in] t0 Init time to apply the pressure
  @param[in] tf Final time to apply the pressure
*/
classHertzianRamp::classHertzianRamp(string type, vector<classElements *> &listOfElements, double P, double T, double mu, double a, vector<double> center, double t0,
                                     double tf) : classNeumannBCs(type) {
    _listOfElements = listOfElements;
    _tf = tf;
    _t0 = t0;
    _a = a;
    _mu = mu;
    _P0 = P*1.5/(3.14*_a*_a);
    _Q0 = _P0 * _mu;
    _c = a * pow( 1 - T/(P*_mu),0.3333);
    _center = center;

};


/*! Constructor for Neumann Pressure BCs Inst
  @param[in] type Type of Neumann BCs
  @param[in] listOfElements Nodes at which the forces are applied
  @param[in] P Pressure
  @param[in] t0 Init time to apply the pressure
  @param[in] tf Final time to apply the pressure
*/
classPressureInst::classPressureInst(string type, vector<classElements *> &listOfElements, double P, double t0,
                                     double tf) : classNeumannBCs(type) {
    _P = P;
    _listOfElements = listOfElements;
    _tf = tf;
    _t0 = t0;

};

/*! Constructor for Neumann Pressure BCs Ramp
  @param[in] type Type of Neumann BCs
  @param[in] listOfElements Nodes at which the forces are applied
  @param[in] P TOTAL normal force
  @param[in] Q TOTAL tangential force
  @param[in] mu Mean value of friction coefficient
  @param[in] muvar Standard deviation of friction coefficient
  @param[in] a Radius of contact
  @param[in] center Coordinates of the center
  @param[in] t0 Init time to apply the pressure
  @param[in] tf Final time to apply the pressure
  @param[in] tf Final time to apply the pressure
  @param[in] nstoch Stochastic dimension
  @param[in] C Third order tensor that computes product between stochastic quantities
  @param[in] flagHaar Indicates whether Haar expansion is used or not
  @param[in] resolution Order of the resolution for Haar expansions
*/
classStochasticHertzianRamp::classStochasticHertzianRamp(string type, vector<classElements *> &listOfElements, double P, double T, double mu, double mu_var, double a, vector<double> center, double t0,
                                     double tf, int nstoch, vector<double> C, bool flagHaar, int resolution, string distribution) : classNeumannBCs(type) {
    _nstoch = nstoch;
    _C = C;
    vector<double> _mu(nstoch,0);
    vector<double> ones(nstoch,0);
    vector<double> mubis(nstoch,0);
    vector<double> c(nstoch,0);
    vector<double> c_final(nstoch,0);
    ones[0] = 1;
    _listOfElements = listOfElements;
    _tf = tf;
    _t0 = t0;
    _a = a;
    _mu[0] = mu;
    _mu[1] = mu_var;
    if(flagHaar){
        vector<double> projection = NormalToHaar(resolution, distribution);
        _mu[1] = 0;
        for(int i=0; i<nstoch; i++){
            _mu[i] = _mu[i] + projection[i]*mu_var;
        }
    }

    mubis = _mu;

    _P0 = P*1.5/(3.14*_a*_a);
    _Q0 = _mu;
    linearCombinationVariableStochastic(_Q0,_P0,0);
    _center = center;
    linearCombinationVariableStochastic(mubis,P/T,0);
    DivisionRandomVariableRandomVariableStochastic(ones, mubis, _C,
                                                   c);
    linearCombinationVariableStochastic(c,-1,1);
    //for(int j=0; j<c.size(); j++){
    //    INFO("La valeur de c[%d] est %g", j, c[j]);
    //}
    cubicRootIntegral(c, c_final,_C);
    linearCombinationVariableStochastic(c_final,a,0);
/*    DL_power_1minusx(c, 0.33, _C);
    linearCombinationVariableStochastic(c,a,0);
    for(int j=0; j<c.size(); j++){
        INFO("La valeur de c[%d] est %g", j, c[j]);
    }*/
    _c = c_final;
};


/*! Constructor for Neumann Current BCs Inst (Extra Dof)
  @param[in] type Type of Neumann BCs
  @param[in] listOfElements Elements at which the current is applied
  @param[in] INeumann Current
  @param[in] t0 Init time to apply the current
  @param[in] tf Final time to apply the current
*/
classCurrentInst::classCurrentInst(string type, vector<classElements *> &listOfElements, double INeumann, double t0,
                                   double tf, int field) : classNeumannBCs(type) {
    _INeumann = INeumann;
    _listOfElements = listOfElements;
    _tf = tf;
    _t0 = t0;
    _fieldIndex = field;
};

/*! Constructor for Neumann Volumetric heat flux BCs (Extra Dof)
  @param[in] type Type of Neumann BCs
  @param[in] qmax maximum heat flux (volumetric)
  @param[in] xc, yc, zc source center
  @param[in] a, b, c maximal axis
  @param[in] t0 Init time to apply the heat flux
  @param[in] tf Final time to apply the heat flux
  @param[in] field extra fielf index
*/
classVolHeatFluxGaussian::classVolHeatFluxGaussian(string type, double qmax, double x0, double y0, double z0,
                                            double vx, double vy, double vz,   
                                            double a, double b, double c,
                                            double t0, double tf, int field) 
                                            : classNeumannBCs(type), _qmax(qmax), _x0(x0), _y0(y0), _z0(z0), 
                                            _vx(vx), _vy(vy), _vz(vz),
                                            _a(a), _b(b), _c(c),
                                            _t0(t0), _tf(tf), _fieldIndex(field) 
{
    INFO("allocate heat source [qmax xc yc zc vx vy vz a b c t0 tf field] =[%g %g %g %g %g %g %g %g %g %g %g %g %d]", qmax,x0,y0,z0,vx,vy,vz,a,b,c,t0,tf,field);

};


/*! Constructor for Neumann Volumetric heat flux BCs (Extra Dof)
  @param[in] type Type of Neumann BCs
  @param[in] listOfElements Elements at which the current is applied
  @param[in] r Heat flux (volumetric)
  @param[in] t0 Init time to apply the heat flux
  @param[in] tf Final time to apply the heat flux
*/
classVolHeatFlux::classVolHeatFlux(string type, vector<classElements *> &listOfElements, double r, double t0,
                                   double tf, int field) : classNeumannBCs(type),_fieldIndex(field) {
    _r = r;
    _listOfElements = listOfElements;
    _tf = tf;
    _t0 = t0;
};

/*! Constructor for Neumann Surface heat flux BCs (Extra Dof)
  @param[in] type Type of Neumann BCs
  @param[in] listOfElements Elements at which the current is applied
  @param[in] HeatNeumann Heat flux (surface)
  @param[in] t0 Init time to apply the heat flux
  @param[in] tf Final time to apply the heat flux
*/
classSurfaceHeatFluxInst::classSurfaceHeatFluxInst(string type, vector<classElements *> &listOfElements, vector<double> HeatNeumann,
                                     double t0, double tf, int field) : classNeumannBCs(type) {
    _HeatNeumann = HeatNeumann;
    _listOfElements = listOfElements;
    _tf = tf;
    _t0 = t0;
    _fieldIndex = field;
};

/*! Constructor for Convection BCs (Extra Dof)
  @param[in] type Type of Neumann BCs
  @param[in] listOfElements Elements at which the current is applied
  @param[in] filmCoefficient Film coefficient or  Convection coefficient
  @param[in] sinkTemp Sink temperature
  @param[in] t0 Init time to apply the heat flux
  @param[in] tf Final time to apply the heat flux
*/
classConvection::classConvection(string type, vector<classElements *> &listOfElements, double filmCoefficient,
                                 double sinkTemp,
                                 double t0, double tf, int field) : classNeumannBCs(type) {

    _convectiveCoef = filmCoefficient;
    _sinkTemp = sinkTemp;
    _listOfElements = listOfElements;
    _tf = tf;
    _t0 = t0;
    _fieldIndex = field;
};

/*! Constructor for Radiation BCs (Extra Dof)
  @param[in] type Type of Neumann BCs
  @param[in] listOfElements Elements at which the current is applied
  @param[in] radiationConstant Radiation constant (emissivity times the Stefan-Boltzmann constant)
  @param[in] sinkTemp Sink temperature (K)
  @param[in] t0 Init time to apply the heat flux
  @param[in] tf Final time to apply the heat flux
*/
classRadiation::classRadiation(string type, vector<classElements *> &listOfElements, double radiationConstant,
                               double sinkTemp, double t0, double tf, int field) : classNeumannBCs(type) {

    _radiationConstant = radiationConstant;
    _sinkTemp = sinkTemp;
    _listOfElements = listOfElements;
    _tf = tf;
    _t0 = t0;
    _fieldIndex = field;
};

classGravity::classGravity(string type, double g) : classNeumannBCs(type) {
    _g = g;
};

classFSI::classFSI(string type, vector<int> &boundary_nodes, vector<vector <double> > boundary_coordinates) : classNeumannBCs(type),
  _boundaryNodes(boundary_nodes), _boundaryCoordinates(boundary_coordinates)
{
  
};

