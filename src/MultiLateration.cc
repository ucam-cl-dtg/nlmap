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

#include "NLMAP/MultiLateration.hh"
#include "NLMAP/LaterationSorter.hh"
#include "NLMAP/NLMAPExceptions.hh"

//-----------------------------------------
// Wrapper for multilateration calculations
//-----------------------------------------
MultiLateration::MultiLateration(REAL *x, 
				 REAL *y, 
				 REAL *z,
				 REAL *d,
				 REAL *sigma,
				 int n) : mLatSort(0){
  for (int i=0; i<n; i++) {
    FitData::Datum datum;
    datum.push_back(x[i]);
    datum.push_back(y[i]);
    datum.push_back(z[i]);
    mLatData.AddDatum(datum,d[i],sigma[i]);
  }
  mLatFunc = new LaterationFunction();
  mLatSort = new LaterationSorter(mLatFunc, &mLatData);
}

MultiLateration::~MultiLateration() {
  if (mLatFunc) delete mLatFunc;
  if (mLatSort) delete mLatSort;
}


//----------------------------------------
// This function takes the lateration data
// and tries to get a sufficiently reliable
// position by repeatedly forming a non-linear
// data model and discarding outliers
//----------------------------------------
XYZData MultiLateration::GetPosition(const int  max_it,
				     const REAL min_delta,
				     const REAL convergence)
{

  mModeler.Model(mLatFunc,
		 &mLatData,
		 mLatSort,
		 max_it,
		 min_delta,
		 convergence);
  XYZData pd;
  pd.x = mLatFunc->GetParams()[0];
  pd.y = mLatFunc->GetParams()[1];
  pd.z = mLatFunc->GetParams()[2];
  pd.sigma = mLatFunc->GetError();
  return pd;
}

