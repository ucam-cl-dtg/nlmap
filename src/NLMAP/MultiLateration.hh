/*
  $Header$
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

#ifndef MULTILATERATION_HH
#define MULTILATERATION_HH

#include <NLMAP/Parameters.hh>
#include <NLMAP/IterativeModeler.hh>
#include <NLMAP/LaterationFunction.hh>
#include <NLMAP/LaterationSorter.hh>
#include <exception>

///
/// Struct to hold position
///
struct XYZData {
  REAL x;
  REAL y;
  REAL z;
  REAL sigma;
};


// Wrapper class for positioning
class MultiLateration {
public:
  ///
  /// Constructor
  /// @param x Array of x co-ordinates
  /// @param y Array of y co-ordinates
  /// @param z Array of z co-ordinates
  /// @param d Array of lateration measures
  /// @param sigma Array of 1-sigma errors on d
  /// @param n Size of input arrays
  ///
  MultiLateration(REAL *x, 
		  REAL *y, 
		  REAL *z,
		  REAL *d,
		  REAL *sigma,
		  int n);
  
  virtual ~MultiLateration();

  ///
  /// Wrapper to use iterative modeler
  /// to calculate a position
  /// @param max_it Maximum number of iterations in each NLM
  /// @param min_delta Minimum delta for NLM to stop
  /// @param convergence Accuracy sought
  /// @return The position
  ///
  XYZData GetPosition(
		      const int  max_it,
		      const REAL min_delta,
		      const REAL convergence);  

  ///
  /// Access to the LaterationFunction
  /// @return Pointer to the LaterationFunction
  ///
  LaterationFunction * GetLaterationFunction() { return mLatFunc; }
  
protected:
  IterativeModeler    mModeler;
  LaterationFunction *mLatFunc;
  LaterationData      mLatData;
  LaterationSorter   *mLatSort;
};

#endif





