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

#ifndef MULTIANGULATION_HH
#define MULTIANGULATION_HH

#include <NLMAP/Parameters.hh>
#include <NLMAP/IterativeModeler.hh>
#include <NLMAP/AngulationFunction.hh>
#include <NLMAP/MagnitudeSorter.hh>
#include <exception>

///
/// Struct to hold result
///
struct XYData {
  REAL x;
  REAL y;
  REAL sigma;
};


///
/// Wrapper class for 2D positioning
/// using multiangulation
///
class MultiAngulation {
public:
  ///
  /// Constructor
  /// @param x Array of x co-ordinates
  /// @param y Array of y co-ordinates
  /// @param b Array of bearing measurements
  /// @param sigma Array of 1-sigma estimates for the b array
  /// @param n Size of the arrays
  ///
  MultiAngulation(REAL *x, 
		  REAL *y, 
		  REAL *b,
		  REAL *sigma,
		  int n);
  
  virtual ~MultiAngulation();

  ///
  /// Wrapper to call the iterative modeler
  /// and extract the result
  /// @param max_it Maximum number of iterations in each NLM
  /// @param min_delta Minimum delta for NLM to stop
  /// @param convergence Accuracy sought
  /// @return The position
  ///
  XYData GetPosition(
		     const int  max_it,
		     const REAL min_delta,
		     const REAL convergence); 

  ///
  /// Access to the AngulationFunction
  /// @return Pointer to the AngulationFunction
  ///
  AngulationFunction * GetAngulationFunction() { return mAngFunc; }
protected:
  IterativeModeler    mModeler;
  AngulationData      mAngData;
  AngulationFunction *mAngFunc;
  MagnitudeSorter    *mAngSort;
};

#endif





