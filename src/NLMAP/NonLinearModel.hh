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

#ifndef NON_LINEAR_MODEL_HH
#define NON_LINEAR_MODEL_HH

#include <NLMAP/Parameters.hh>
#include <exception>
#include <vector>

/**
 * Class to handle input data
 */
class FitData {
 public:
  FitData(const int numFields);
  virtual ~FitData();

  typedef std::vector<REAL> Datum;

  // Get number of valid measurments in set
  virtual int   GetInputDataSize();
  // Get measurement for index idx
  virtual REAL  GetMeasurement(const int idx);
  // Get the error estimate of measurement idx
  virtual REAL  GetSigma(const int idx);
  // Get validity of specific measurement
  virtual bool  GetValidity(const int idx);
  // Get the std error of the model parameters
  virtual void  Invalidate(const int n);
  // Get the raw data
  virtual REAL  GetData(const int field, const int idx);
  // Add some new data
  virtual void  AddDatum (Datum &data,
			  REAL measurement,
			  REAL error);

protected:
  // Convert an index into valid data
  // to an index into all data (inc. invalid)
  int  ActualIndex(int j);

  int                  mNumFields;
  std::vector< Datum > mData;
  std::vector< REAL >  mMeasurements;
  std::vector< REAL >  mSigma;
  std::vector< bool >  mValidity;
};



/**
 * Pure virtual class overridden
 * with a custom function and parameters
 */
class FitFunction {
 public:
  FitFunction(const int nparams, FitData *fd);
  virtual ~FitFunction();
  // Return the residual between two values
  virtual REAL CalculateResidual(REAL val1, REAL val2);// {return val1-val2;}
  // Return the number of parameters in the fitting function
  virtual int   GetNumParameters();
  // Initialise the parameters
  virtual void  InitialiseParameters()=0;
  // Return the parameter estimates
  virtual REAL* GetParams();
  // Get derivative wrt parameters for index idx
  virtual REAL  GetCurrentDerivative(const int idx);
  // Evaluate the fitting function
  virtual REAL  Evaluate(const int idx, REAL parameters[])=0;
  // Get the std error of the model parameters
  virtual REAL  GetError();
  // Set the std error of the model parameters;
  virtual void  SetError(const REAL s);
protected:
  int      mNumParams;
  REAL     mError;
  REAL    *mParams;
  REAL    *mDeriv;
  FitData *mData;
};









/**
 * Pure virtual class to override with
 * your own function 
 */
//  class FittingFunction {
//    public:
//    virtual int   GetInputDataSize()=0;
//    // Get measurement for index idx
//    virtual REAL  GetMeasurement(const int idx)=0;
//    // Get the error estimate of measurement idx
//    virtual REAL  GetSigma(const int idx)=0;
//    // Get validity of specific measurement
//    virtual bool  GetValidity(const int idx)=0;
//    // Get the std error of the model parameters
//    virtual REAL  GetError()=0;
//    // Set the std error of the model parameters;
//    virtual void  SetError(const REAL s)=0;
//    // Return the residual between two values
//    virtual REAL CalculateResidual(REAL val1, REAL val2)=0;
//    // Return the number of parameters in the fitting function
//    virtual int   GetNumParameters()=0;
//    // Initialise the parameters
//    virtual void  InitialiseParameters()=0;
//    // Return the parameter estimates
//    virtual REAL* GetParams()=0;
//    // Get derivative wrt parameters for index idx
//    virtual REAL  GetCurrentDerivative(const int idx)=0;
//    // Evaluate the fitting function
//    virtual REAL  Evaluate(const int x, REAL parameters[])=0;
//    // Invalidate a datum
//    virtual void  Invalidate(int n)=0;
//  };


class NonLinearModel {
public:
  NonLinearModel(FitFunction *ff, FitData *fd);
  virtual ~NonLinearModel();

  REAL GetStdErr();


  void Fit(int max_it,
	   REAL min_delta);

protected:

  // Perform a single iteration
  int SingleMarquardtIteration();

  // Find alpha, beta, chisq
  void ComputeSupportData(REAL *params);

  FitFunction    *mFitFunction;
  FitData        *mFitData;
  int             mDataSize;
  int             mNumParams;
  REAL            mChiSq;
  REAL            mLambda;
  REAL          **mAlpha;
  REAL          **mAlphadash;
  REAL           *mBeta;
  REAL           *mParams;
  REAL          **mDa;
  REAL          **mSingleDa;
  bool            mLastCall;
};

#endif


