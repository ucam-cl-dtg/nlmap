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

#ifndef ANGULATION_FUNCTION_HH
#define ANGULATION_FUNCTION_HH

#include <NLMAP/NonLinearModel.hh>
#include <NLMAP/NLMAPExceptions.hh>

#include <vector>


///
/// A specific data store for 2D 
/// positioning by angulation
///
class AngulationData : public FitData 
{
public:

  ///
  /// Constructor
  ///
  AngulationData();

  ///
  /// Destructor
  ///
  virtual ~AngulationData();

  ///
  /// Get the X co-ordinate of base
  /// station
  /// @param i The station index
  ///
  REAL GetX(const int i);

  ///
  /// Get the Y co-ordinate of base
  /// station
  /// @param i The station index
  ///
  REAL GetY(const int i);
};



///
/// Specification of the parameters
/// and the fit function for 2D
/// positioning by angulation
///
class AngulationFunction : public FitFunction
{
public:

  ///
  /// Constructor
  ///
  AngulationFunction(const int nparams, AngulationData *d);
  
  ///
  /// Destructor
  ///
  virtual ~AngulationFunction();

  ///
  /// Initialise the parameters
  ///
  virtual void  InitialiseParameters();
  
  ///
  /// Evaluate the fitting function
  /// for a given datum
  /// @param idx Index of the datum
  /// @param parameters Array of parameter
  /// values to evaluate against
  ///
  virtual REAL  Evaluate(const int idx, REAL parameters[]);
  
  ///
  /// Special residual calc needed
  /// to account for modulo-2pi nature
  /// of bearings.
  /// @param b First bearing
  /// @param b2 Second bearing
  ///
  virtual REAL CalculateResidual(REAL b, REAL b2);
};



#endif//ANGULATION_FUNCTION_HH

