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

#ifndef LATERATION_FUNCTION_HH
#define LATERATION_FUNCTION_HH

#include <NLMAP/NonLinearModel.hh>
#include <NLMAP/NLMAPExceptions.hh>
#include <vector>



///
/// A data store for lateration results
/// in 3D
///
class LaterationData : public FitData {
public:
  ///
  /// Constructor
  ///
  LaterationData();

  ///
  /// Destructor
  ///
  virtual ~LaterationData();
  
  ///
  /// Get X co-ordinate for
  /// the i'th valid datum
  /// @param i Index of valid data
  ///
  REAL GetX(const int i);

  ///
  /// Get Y co-ordinate for
  /// the i'th valid datum
  /// @param i Index of valid data
  ///
  REAL GetY(const int i);
  
  ///
  /// Get Z co-ordinate for
  /// the i'th valid datum
  /// @param i Index of valid data
  ///
  REAL GetZ(const int i);
};




///
/// This is a specific fit function
/// for lateration data. It allows
/// an (x,y,z) estimate given a series
/// of ranging measures from beacons
/// at known positions.
///
class LaterationFunction : public FitFunction {
public:

  ///
  /// Constructor
  /// @param nparams Number of parameters (3?)
  /// @param d Pointer to the lateration data
  ///
  LaterationFunction(const int nparams, LaterationData *d);

  virtual ~LaterationFunction();

  ///
  /// Initialise the parameters by
  /// averaging the input locations
  /// to make an initial seed guess
  ///
  virtual void  InitialiseParameters();

  ///
  /// Evaluate the fitting function
  ///
  virtual REAL  Evaluate(const int idx, REAL parameters[]);

};



#endif


