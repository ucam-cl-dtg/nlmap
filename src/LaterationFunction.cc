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


LaterationFunction::LaterationFunction(const int nparams, LaterationData *d)
  : FitFunction(nparams,d) {}

LaterationFunction::~LaterationFunction() {}


//----------------------------------------
// Return the number of parameters in the fitting function
//----------------------------------------
int LaterationFunction::GetNumParameters() {return 3;}



//----------------------------------------
// Initialise the parameters ready for 
// regression
//----------------------------------------
void LaterationFunction::InitialiseParameters()
{
  LaterationData *data = dynamic_cast<LaterationData *>(mData);
  // Guess at the average position
  for (int i=0; i<3; i++) mParams[i]=0.0;
  REAL z_min = data->GetZ(0);
  REAL d_min = data->GetMeasurement(0);
  for (int i=0; i<data->GetInputDataSize(); i++) {
    mParams[0] += data->GetX(i);
    mParams[1] += data->GetY(i);
    if (data->GetMeasurement(i) < d_min) d_min=data->GetMeasurement(i);
    if (data->GetZ(i) < z_min) z_min=data->GetZ(i);
  }
  mParams[0] /= data->GetInputDataSize();
  mParams[1] /= data->GetInputDataSize();
  mParams[2]  = z_min-d_min;
  
  mParams[0] = 2.0/3.0+0.5;
  mParams[1] = 4.0/3.0-0.5;
  mParams[2]  = 1.0/3.0+0.5;
}


//----------------------------------------
// Evaluate the current fit function
//----------------------------------------
REAL LaterationFunction::Evaluate(const int i, REAL *parameters) 
{
  
  LaterationData *data = dynamic_cast<LaterationData *>(mData);
  REAL d = sqrt(
		(data->GetX(i)-parameters[0])*(data->GetX(i)-parameters[0]) +
		(data->GetY(i)-parameters[1])*(data->GetY(i)-parameters[1]) +
		(data->GetZ(i)-parameters[2])*(data->GetZ(i)-parameters[2])
		);
  mDeriv[0] = (parameters[0]-data->GetX(i))/d;
  mDeriv[1] = (parameters[1]-data->GetY(i))/d;
  mDeriv[2] = (parameters[2]-data->GetZ(i))/d;
  return d;
}


