//
//
// File authors: see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file timeSteping.cpp
  \brief This file contains all functions related the time stepping procedure
*/

#include "timeStepping.h"
#include "solvers.h"
#include "commandLine.h"

AdaptiveTimeStep* AdaptiveTimeStep::createAdaptiveTimeStep(const char what[], const std::vector<double>& data)
{
    if (strcmp(what,"POWER") == 0 )
    {
        if (data.size() <2 )
        {
            ERROR("The number of parameters must be equal to 2 with PowerLawAdaptiveTimeStep");
            exit(-1);
        }
        return new PowerLawAdaptiveTimeStep(std::round(data[0]),data[1]);
    }
    else
    {
        ERROR("AdaptiveTimeStep is not defined with %s",what);
        exit(-1);
    }
};

PowerLawAdaptiveTimeStep::PowerLawAdaptiveTimeStep(int nIterOpt, double n): AdaptiveTimeStep(), _nIterOpt(nIterOpt), _power(n)
{
    if (GeneralOptions::commRank == 0)
        INFO("Adaptive time step with a power law, number optimal iteration = %d, power = %g",_nIterOpt,_power);
}

double PowerLawAdaptiveTimeStep::getTimeStep(int Nconverged, double dtPrev) const
{
    double fact = (double)_nIterOpt/(double)Nconverged;
    double dt = dtPrev*pow(fact,_power);
    if (GeneralOptions::commRank == 0)
    {
        if (fact >1)
        {
            INFO("The time step increase from %g to %g",dtPrev, dt);
        }
        else if (fact <1)
        {
            INFO("The time step decrease from %g to %g",dtPrev, dt);
        }
    }
    return dt;
}


TimeStepping::TimeStepping(double tol):
    _startTime(0.),_endTime(0.), _numSteps(0),_tol(tol),
    _lastTime(0.),_timeStep(0.), _timeStepMax(1e10),
    _times(1,0),_nsteps(),_forced(false){}

TimeStepping::TimeStepping(const TimeStepping& src):
    _startTime(src._startTime),
    _endTime(src._endTime),
    _numSteps(src._numSteps),
    _tol(src._tol),
    _times(src._times),
    _nsteps(src._nsteps),
    _timeStep(src._timeStep),
    _timeStepMax(src._timeStepMax),
    _lastTime(src._lastTime),
    _forced(src._forced)
{
}

TimeStepping& TimeStepping::operator = (const TimeStepping& src)
{
    _startTime = src._startTime;
    _endTime = src._endTime;;
    _numSteps = src._numSteps;
    _tol = src._tol;
        //
    _times = src._times;
    _nsteps = src._nsteps;
        
         //
    _timeStep = src._timeStep;
    _timeStepMax = src._timeStepMax;
    _lastTime = src._lastTime;
    
    //
    _forced = src._forced;
        
    return *this;
}

    
TimeStepping::~TimeStepping(){}

void TimeStepping::reset()
{
    _startTime = 0;
    _endTime = 0.;
    _numSteps = 0;
    _tol = 1e-6;
    _lastTime = 0.;
    _timeStep = 0.;
    _times.clear();
    _times.resize(1,0);
    _nsteps.clear();
    _forced = false;
}

void TimeStepping::setMaximumTimeStep(double dtMax)
{
    _timeStepMax = dtMax;
}

void TimeStepping::setTime(double lastTime, double dt)
{
    _lastTime = lastTime;
    _timeStep = dt;
}

void TimeStepping::setTimeStep(double dt)
{
    _timeStep = dt;
    if (_forced) return;
    
    _timeStep = std::min(dt,_timeStepMax);
    
    // update time step
    for (int i=0; i<_nsteps.size(); i++)
    {
        double lowerBound = _times[i];
        double upperBound = _times[i+1];
        
        if (fabs(_lastTime-lowerBound) < _tol*_endTime)
        {
            // start in interval
            _timeStep = std::min(_timeStepMax,(upperBound-lowerBound)/_nsteps[i]);
            break;
        }
        else if (_lastTime < upperBound-_tol*_endTime and _lastTime+_timeStep > upperBound-_tol*_endTime)
        {
            // truncate
            _timeStep = upperBound - _lastTime;
            break;
        };
    } 
    // truncate
    if (_lastTime + _timeStep > _endTime-_tol*_endTime)
    {
        _timeStep = _endTime - _lastTime;
    }
}
void TimeStepping::forceSteppingScheme(bool fl)
{
    _forced = fl;
}

void TimeStepping::setStartTime(double st)
{
    _startTime = st;
    _times[0] = _startTime;
}

void TimeStepping::setEndTime(double endt)
{
    _endTime = endt;
}
        
void TimeStepping::setNumberOfSteps(int nstep)
{
    _numSteps = nstep;
}
void TimeStepping::setNumberOfSteps(double time, int nstep)
{
    _times.push_back(time);
    _nsteps.push_back(nstep);
}

void TimeStepping::initFromOptions()
{
    if (_nsteps.size() >0)
    {
            if (fabs(_endTime - _times.back()) >_tol)
            {
                ERROR("time stepping plan is not correctly defined, last time must be equal to end time");
                exit(-1);
            }   
            _numSteps = 0;
            for (int i=0; i< _nsteps.size(); i++)
            {
                _numSteps += _nsteps[i];
            }
    }
    else
    {
        _times.push_back(_endTime);
        _nsteps.push_back(_numSteps);
    }
    
    // initial
    _lastTime = _startTime;
    if (_numSteps >0)
    {
      _timeStep = std::min(_timeStepMax,(_times[1] - _times[0])/_nsteps[0]);
    }
    else
    {
      _timeStep = 0.;
    }
}

void TimeStepping::nextStep()
{
    // update last time
    _lastTime += _timeStep;
    
    if (_forced) return;
    
    _timeStep = std::min(_timeStep,_timeStepMax);
    // update time step
    for (int i=0; i<_nsteps.size(); i++)
    {
        double lowerBound = _times[i];
        double upperBound = _times[i+1];
        
        if (fabs(_lastTime-lowerBound) < _tol*_endTime)
        {
            // start in interval
            _timeStep = std::min(_timeStepMax,(upperBound-lowerBound)/_nsteps[i]);
            break;
        }
        else if (_lastTime < upperBound-_tol*_endTime and _lastTime+_timeStep > upperBound-_tol*_endTime)
        {
            // truncate
            _timeStep = upperBound - _lastTime;
            break;
        };
    } 
    // truncate
    if (_lastTime + _timeStep > _endTime-_tol*_endTime)
    {
        _timeStep = _endTime - _lastTime;
    }
}

void TimeStepping::modifyTimeStepByFactor(double factor)
{
    _timeStep *= factor;
}

bool TimeStepping::endCheck() const
{
    if (fabs(_lastTime - _endTime) < _tol) 
    {
        return false;
    }
    else
    {
        return true;
    }
 
};