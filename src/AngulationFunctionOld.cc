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

#include "NLMAP/AngulationFunction.hh"
#include <cmath>

AngulationFunction::AngulationFunction() : 
  mInit(false), mStdErr(-1.0), mValidity(0) {};

AngulationFunction::~AngulationFunction() {
  if (mValidity) delete mValidity;
}

void AngulationFunction::AddDatum(REAL x, REAL y,
				  REAL bearing,
				  REAL sigma)
{
  if (mInit) throw InitialisationError();
  if (sigma<=0.0 || bearing<0.0 || bearing>2.0*M_PI) throw InvalidData();
  mRx.push_back(x);
  mRy.push_back(y);
  mBearing.push_back(bearing);
  mSigma.push_back(sigma);
}

//----------------------------------------
// Get the number of VALID measurements
// left in the dataset
//----------------------------------------
int AngulationFunction::GetInputDataSize() {
  if (!mInit) throw Uninitialised();
  int n=0;
  for (int i=0; i<mBearing.size(); i++) {
    n+=mValidity[i];
  }
  return n;
}

//----------------------------------------
// Initialise the parameters ready for 
// regression
//----------------------------------------
void AngulationFunction::InitialiseParameters()
{

  if (!mValidity) {
    mValidity = new bool[mBearing.size()];
    for (int i=0;i<mBearing.size(); i++) mValidity[i]=true;
    mInit=true;
  }

  // Use the first two measurements to get a start position
  REAL diff1y = cos(GetMeasurement(0));
  REAL diff1x = (GetMeasurement(0)<M_PI) ? 
    sqrt(1-diff1y*diff1y):-sqrt(1-diff1y*diff1y);
  REAL diff2y = cos(GetMeasurement(1));
  REAL diff2x = (GetMeasurement(1)<M_PI) ? 
    sqrt(1-diff2y*diff2y):-sqrt(1-diff2y*diff2y);
  
  REAL mu = ((GetX(0)-GetX(1))*diff1y - (GetY(0)-GetY(1))*diff1x) /
    (diff2x*diff1y-diff2y*diff1x);
  
  mParams[0]=GetX(1)+mu*diff2x;
  mParams[1]=GetY(1)+mu*diff2y;
}

//----------------------------------------
// Translate an index for a VALID 
// measurement into an index into the full
// dataset
//----------------------------------------
int AngulationFunction::ActualIndex(int i) {
  // return i;
  int n=0;
  for (int j=0; j<mBearing.size(); j++) {
    if (n==i && mValidity[j]) return j;
    if (mValidity[j])n++;
  }
  return -1;
}


REAL AngulationFunction::Evaluate(const int i, REAL *parameters) 
{
  if (!mInit) throw Uninitialised();
  int idx = ActualIndex(i);
  
  // First calculate the bearing of the estimate
  // Cross product of (0,1,0) with (px-rix,py-riy,0)
  // gives
  REAL pminusrix = parameters[0]-mRx[idx];
  REAL pminusriy = parameters[1]-mRy[idx];
  REAL pminusri = sqrt(pminusrix*pminusrix+pminusriy*pminusriy);

  REAL b = acos(pminusriy/pminusri);
  if (pminusrix<0.0) b=2.0*M_PI-b;

  // Now the derivatives
  mDeriv[1]=-pminusrix/pminusri;
  mDeriv[0]=pminusriy/pminusri;

  return b;
}

// Specialised residual function because everything
// is modulo 2PI
REAL AngulationFunction::CalculateResidual(REAL b, REAL b2)
{
  REAL res1 = b-b2;
  REAL res2;
  if (b>b2) res2=(b+2.0*M_PI)-b2;
  else res2=b-(b2+2.0*M_PI);
  
  REAL res = abs(res1)<=abs(res2)?abs(res1):abs(res2);
  return res;
}
