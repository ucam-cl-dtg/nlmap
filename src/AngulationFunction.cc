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

AngulationData::AngulationData()
  : FitData(2){}

AngulationData::~AngulationData() {}


REAL AngulationData::GetX(const int i) {
  return (mData[ActualIndex(i)])[0];
}


REAL AngulationData::GetY(const int i)
{
  return (mData[ActualIndex(i)])[1];
}


AngulationFunction::AngulationFunction()
  :  FitFunction(2) {};


AngulationFunction::~AngulationFunction() {}


//----------------------------------------
// Initialise the parameters ready for 
// regression
//----------------------------------------
void AngulationFunction::InitialiseParameters(FitData *fd)
{

  AngulationData *data = dynamic_cast<AngulationData *>(fd);

  // Use the first two measurements to get a start position
  REAL diff1y = cos(data->GetMeasurement(0));
  REAL diff1x = (data->GetMeasurement(0)<M_PI) ? 
    sqrt(1-diff1y*diff1y):-sqrt(1-diff1y*diff1y);
  REAL diff2y = cos(data->GetMeasurement(1));
  REAL diff2x = (data->GetMeasurement(1)<M_PI) ? 
    sqrt(1-diff2y*diff2y):-sqrt(1-diff2y*diff2y);
  
  REAL mu = ((data->GetX(0)-data->GetX(1))*diff1y - 
	     (data->GetY(0)-data->GetY(1))*diff1x) /
    (diff2x*diff1y-diff2y*diff1x);
  
  mParams[0]=data->GetX(1)+mu*diff2x;
  mParams[1]=data->GetY(1)+mu*diff2y;
}



REAL AngulationFunction::Evaluate(const int i, REAL *parameters, FitData *fd) 
{
  AngulationData *data = dynamic_cast<AngulationData *>(fd);

  // First calculate the bearing of the estimate
  // Cross product of (0,1,0) with (px-rix,py-riy,0)
  // gives
  REAL pminusrix = parameters[0]-data->GetX(i);
  REAL pminusriy = parameters[1]-data->GetY(i);
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
  
  REAL res = fabs(res1)<=fabs(res2)?fabs(res1):fabs(res2);
  return res;
}
