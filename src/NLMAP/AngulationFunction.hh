/*
  Copyright (C) 2004 Robert K. Harle

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Email: robert.harle@achilles.org
  Email: rkh23@cantab.net
*/

#ifndef ANGULATION_FUNCTION_HH
#define ANGULATION_FUNCTION_HH

#include <NLMAP/NonLinearModel.hh>
#include <NLMAP/NLMAPExceptions.hh>

#include <vector>


class AngulationData : public FitData 
{
public:
  AngulationData();
  virtual ~AngulationData();
  REAL GetX(const int i);
  REAL GetY(const int i);
};

class AngulationFunction : public FitFunction
{
public:
  AngulationFunction(const int nparams, AngulationData *d);
  virtual ~AngulationFunction();
  // Initialise the parameters
  virtual void  InitialiseParameters();
  // Evaluate the fitting function
  virtual REAL  Evaluate(const int idx, REAL parameters[]);
  // Special residual calc
  virtual REAL CalculateResidual(REAL b, REAL b2);
};







//  class AngulationFunction : public FitFunction {
//  public:
//    AngulationFunction();
//    virtual ~AngulationFunction();

//    int   GetInputDataSize();
//    int   GetNumParameters() { return 2;}
//    REAL* GetParams() { return mParams; }
//    REAL  GetCurrentDerivative(const int idx) { return mDeriv[idx]; }
//    REAL  GetMeasurement(const int idx) { return mBearing[ActualIndex(idx)];} 
//    REAL  GetSigma(const int idx) { return mSigma[ActualIndex(idx)]; }
//    bool  GetValidity(const int idx) { return mValidity[ActualIndex(idx)]; }
//    REAL  Evaluate(const int idx, REAL *parameters);
//    REAL  GetError() { return mStdErr; }
//    void  SetError(const REAL err) { mStdErr=err; }
//    REAL  CalculateResidual(REAL a, REAL b);

//    REAL  GetX(const int i) { return mRx[ActualIndex(i)]; }
//    REAL  GetY(const int i) { return mRy[ActualIndex(i)]; }

//    void  Invalidate(const int i) {mValidity[ActualIndex(i)]=false;}
//    void  InitialiseParameters();
//    void  AddDatum(REAL x, REAL y,
//  		 REAL bearing,
//  		 REAL sigma);

//  private:
//    int ActualIndex(const int i);

//    std::vector<REAL>    mRx;
//    std::vector<REAL>    mRy;
//    std::vector<REAL>    mSigma;
//    std::vector<REAL>    mBearing;
//    bool                *mValidity;
//    REAL                 mParams[2];
//    REAL                 mDeriv[2];
//    bool                 mInit;
//    REAL                 mStdErr;
//  };

#endif
