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

#ifndef LATERATION_FUNCTION_HH
#define LATERATION_FUNCTION_HH

#include <NLMAP/NonLinearModel.hh>
#include <NLMAP/NLMAPExceptions.hh>
#include <vector>


//--------------------------------------
// This is a specific fit function
// for lateration data. It allows
// an (x,y,z) estimate given a series
// of ranging measures from beacons
// at known positions.
//
// Note that data is never removed
// as a model iterates: only invalidated.
// This requires translation of indexes
// which is done using the private
// ActualIndex() member function
//--------------------------------------
//  class LaterationFunction : public FitFunction {
//  public:
//    LaterationFunction();
//    virtual ~LaterationFunction();

//    // These are the pure virtual functions we must provide
//    int   GetInputDataSize();
//    int   GetNumParameters() { return 3;}
//    REAL* GetParams() { return mParams; }
//    REAL  GetCurrentDerivative(const int idx) { return mDeriv[ActualIndex(idx)]; }
//    REAL  GetMeasurement(const int idx) { return mD[ActualIndex(idx)];} 
//    REAL  GetSigma(const int idx) { return mSigma[ActualIndex(idx)]; }
//    bool  GetValidity(const int idx) { return mValidity[ActualIndex(idx)]; }
//    REAL  Evaluate(const int idx, REAL *parameters);
//    REAL  GetError() { return mStdErr; }
//    void  SetError(const REAL err) { mStdErr=err; }

//      // These are specific functions for laterations
//    void  Invalidate(const int i) {mValidity[ActualIndex(i)]=false;}
//    REAL  GetX(const int i) { return mX[ActualIndex(i)]; }
//    REAL  GetY(const int i) { return mY[ActualIndex(i)]; }
//    REAL  GetZ(const int i) { return mZ[ActualIndex(i)]; }

//    void  InitialiseParameters();
//    void  AddDatum(REAL x, REAL y, REAL z, 
//  		 REAL laterationMeasurement,
//  		 REAL sigma);
  
//  private:
//    int ActualIndex(const int i);

//    std::vector<REAL>    mX;
//    std::vector<REAL>    mY;
//    std::vector<REAL>    mZ;
//    std::vector<REAL>    mD;
//    std::vector<REAL>    mSigma;
//    bool                *mValidity;
//    REAL                 mParams[3];
//    REAL                 mDeriv[3];
//    bool                 mInit;
//    REAL                 mStdErr;
//  };


class LaterationData : public FitData {
public:
  LaterationData();
  virtual ~LaterationData();
  REAL GetX(const int i);
  REAL GetY(const int i);
  REAL GetZ(const int i);
};



class LaterationFunction : public FitFunction {
public:
  LaterationFunction(const int nparams, LaterationData *d);
  virtual ~LaterationFunction();
  // Initialise the parameters
  virtual void  InitialiseParameters();
  // Evaluate the fitting function
  virtual REAL  Evaluate(const int idx, REAL parameters[]);
  virtual int GetNumParameters();

};



#endif


