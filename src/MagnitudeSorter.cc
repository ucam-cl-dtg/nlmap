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

#include "NLMAP/MagnitudeSorter.hh"
#include "NLMAP/NLMAPExceptions.hh"

#include <cmath>


// This residual sorter class
// sorts residuals by magnitude
int MagnitudeSorter::GetMaxResidualIndex() {

  int n = mData->GetInputDataSize();

  int  maxidx=0;
  REAL maxres=0.0;
  // Calculate the residuals
  for (int i=0; i<n; i++) {
    REAL b  = mFunc->Evaluate(i,mFunc->GetParams(),mData);
    REAL b2 = mData->GetMeasurement(i);

    REAL res = mFunc->CalculateResidual(
					   mFunc->Evaluate(i,
							   mFunc->GetParams(),
							   mData
							   ), 
					   mData->GetMeasurement(i)
					   );
    // Want to know the magnitude of the residual only
    res = fabs(res);
    if (res>maxres) {
      maxres=res;
      maxidx=i;
    }

  }
  return maxidx;
}
