//
//
// File authors: see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file timeSteping.h
  \brief This file contains all functions related the time stepping procedure
*/


#ifndef _TIMESTEPPING_H_
#define _TIMESTEPPING_H_

#include <vector>

class AdaptiveTimeStep
{
    public:
        virtual ~AdaptiveTimeStep(){}
        virtual double getTimeStep(int Nconverged, double dtPrev) const =0;
        static AdaptiveTimeStep* createAdaptiveTimeStep(const char what[], const std::vector<double>& data);
    
};

class PowerLawAdaptiveTimeStep : public AdaptiveTimeStep
{
    protected:
        int _nIterOpt;
        double _power;
    
    public:
        PowerLawAdaptiveTimeStep(int nIterOpt, double n);
        virtual ~PowerLawAdaptiveTimeStep(){}
        virtual double getTimeStep(int Nconverged, double dtPrev) const;
};

class TimeStepping
{
    protected:
        double _startTime, _endTime;
        int _numSteps;
        double _tol;
         //
        double _timeStep;
        double _lastTime;
        double _timeStepMax; // maximum time step possible
        //
        std::vector<double> _times;
        std::vector<int> _nsteps;
        //
        bool _forced; // when FSI is employed, the time stepping if forced
        
    public:
        TimeStepping(double tol=1e-12);
        TimeStepping(const TimeStepping& src);
        TimeStepping& operator = (const TimeStepping& src);
        ~TimeStepping();
        
        void reset();
        void setStartTime(double st);
        void setEndTime(double endt);
        void setNumberOfSteps(int nstep);
        void setNumberOfSteps(double time, int nstep);
        
        void setMaximumTimeStep(double dtMax);
        void setTimeStep(double dt);
        void setTime(double lastTime, double dt);
        void forceSteppingScheme(bool fl);
        
        void initFromOptions();
        double getTimeStep() const {return _timeStep;};
        double getLastTime() const {return _lastTime;};
        double getTimeRun() const {return _lastTime+_timeStep;};
        double getNumSteps() const {return _numSteps;};
        double getStartTime() const {return _startTime;};
        double getEndTime() const { return _endTime;};
        void modifyTimeStepByFactor(double factor);
        void nextStep();
        bool endCheck() const;
        
        TimeStepping* clone() const {return new TimeStepping(*this);};
};
#endif // _TIMESTEPPING_H_