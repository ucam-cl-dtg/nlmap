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

#include "NLMAP/MultiAngulation.hh"
#include "NLMAP/MagnitudeSorter.hh"
#include "NLMAP/NLMAPExceptions.hh"

//-----------------------------------------
// Wrapper for multilateration calculations
//-----------------------------------------
MultiAngulation::MultiAngulation(REAL *x, 
				 REAL *y, 
				 REAL *d,
				 REAL *sigma,
				 int n) : mAngSort(0){
  for (int i=0; i<n; i++) {
    FitData::Datum datum;
    datum.push_back(x[i]);
    datum.push_back(y[i]);
    mAngData.AddDatum(datum,d[i],sigma[i]);
  }
  mAngFunc = new AngulationFunction();
  mAngSort = new MagnitudeSorter(mAngFunc, &mAngData);
}

MultiAngulation::~MultiAngulation() {
  if (mAngFunc) delete mAngFunc;
  if (mAngSort) delete mAngSort;
}


//----------------------------------------
// This function takes the lateration data
// and tries to get a sufficiently reliable
// position by repeatedly forming a non-linear
// data model and discarding outliers
//----------------------------------------
XYData MultiAngulation::GetPosition(const int  max_it,
				       const REAL min_delta,
				       const REAL convergence)
{
  mModeler.Model(mAngFunc,
		 &mAngData,
		 mAngSort,
		 max_it,
		 min_delta,
		 convergence);

  XYData pd;
  pd.x = mAngFunc->GetParams()[0];
  pd.y = mAngFunc->GetParams()[1];
  pd.sigma = mAngFunc->GetError();
  return pd;
}

