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

#include "NLMAP/LaterationFunction.hh"

#include <cmath>

using namespace std;

LaterationFunction::LaterationFunction() : 
  mInit(false), mStdErr(-1.0), mValidity(0) {};

LaterationFunction::~LaterationFunction() {
  if (mValidity) delete mValidity;
}



//----------------------------------------
// Add a measurement to the dataset
// Once Initialisation has occured
// this should NOT be called
//----------------------------------------
void LaterationFunction::AddDatum(REAL x, REAL y, REAL z, 
			     REAL laterationMeasurement,
			     REAL sigma)
{
  if (mInit) throw InitialisationError();
  if (sigma<=0.0 || laterationMeasurement<0.0) throw InvalidData();
  mX.push_back(x);
  mY.push_back(y);
  mZ.push_back(z);

  mD.push_back(laterationMeasurement);
  mSigma.push_back(sigma);
}



//----------------------------------------
// Get the number of VALID measurements
// left in the dataset
//----------------------------------------
int LaterationFunction::GetInputDataSize() {
  if (!mInit) throw Uninitialised();
  int n=0;
  for (int i=0; i<mD.size(); i++) {
    n+=mValidity[i];
  }
  return n;
}



//----------------------------------------
// Initialise the parameters ready for 
// regression
//----------------------------------------
void LaterationFunction::InitialiseParameters()
{
  if (!mValidity) {
    mInit=true;

    mValidity = new bool[mX.size()];
    for (int i=0;i<mD.size(); i++) mValidity[i]=true;
  }


  // Guess at the average position
  for (int i=0; i<3; i++) mParams[i]=0.0;
  REAL z_min = GetZ(0);
  REAL d_min = GetMeasurement(0);
  for (int i=0; i<GetInputDataSize(); i++) {
    mParams[0] += GetX(i);
    mParams[1] += GetY(i);
    if (GetMeasurement(i) < d_min) d_min=GetMeasurement(i);
    if (GetZ(i) < z_min) z_min=GetZ(i);
  }
  mParams[0] /= GetInputDataSize();
  mParams[1] /= GetInputDataSize();
  mParams[2]  = z_min-d_min;
}



//----------------------------------------
// Translate an index for a VALID 
// measurement into an index into the full
// dataset
//----------------------------------------
int LaterationFunction::ActualIndex(int i) {
  // return i;
  int n=0;
  for (int j=0; j<mD.size(); j++) {
    if (n==i && mValidity[j]) return j;
    if (mValidity[j])n++;
  }
  return -1;
}



//----------------------------------------
// Evaluate the current fit function
//----------------------------------------
REAL LaterationFunction::Evaluate(const int i, REAL *parameters) 
{
  if (!mInit) throw Uninitialised();
  int idx = ActualIndex(i);
  REAL d = sqrt(
		(mX[idx]-parameters[0])*(mX[idx]-parameters[0]) +
		(mY[idx]-parameters[1])*(mY[idx]-parameters[1]) +
		(mZ[idx]-parameters[2])*(mZ[idx]-parameters[2])
		);
  mDeriv[0] = (parameters[0]-mX[idx])/d;
  mDeriv[1] = (parameters[1]-mY[idx])/d;
  mDeriv[2] = (parameters[2]-mZ[idx])/d;
  return d;
}


