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

#ifndef ITERATIVE_MODELER_HH
#define ITERATIVE_MODELER_HH

#include <NLMAP/Parameters.hh>
#include <NLMAP/NonLinearModel.hh>
#include <NLMAP/NLMAPExceptions.hh>

///
/// Pure virtual class to override with your
/// specifics for ordering the model residuals
/// and discarding outliers
///
class ResidualSorter {
public:
  ///
  /// Constructor
  /// @param ff Pointer to the fitting function
  /// @param fd Pointer to the fitting data
  ///
  ResidualSorter(FitFunction *ff, FitData *fd);

  ///
  /// Destructor
  ///
  virtual ~ResidualSorter();

  ///
  /// Calculate and order the residuals
  /// using some metric and return the index
  /// of the datum to discard
  /// @return Index of maximum residual datum
  ///
  virtual int GetMaxResidualIndex()=0;
protected:
  FitFunction  *mFunc;
  FitData      *mData;
};



///
/// Wrapper class for iteratively modelling
/// i.e. forming NLM of data, evaluating
/// the error, and discarding the outliers
/// and repeating until either sufficient
/// accuracy is achieved or there is not enough
/// data left to form a NLM
///
class IterativeModeler {
public:

  ///
  /// Empty constructor
  ///
  IterativeModeler() {};


  ///
  /// Initiate an iterative model
  /// of some data using a specific
  /// function and sorting residuals
  /// according to some sorter
  /// @param ff Pointer to FitFunction object
  /// @param fd Pointer to FitData object
  /// @param rs Pointer to ResidualSorter object
  /// @param max_iterations Maximum iterations allowed FOR EACH NLM
  /// @param min_delta Minimum delta to achieve with each NLM
  /// @param convergence The accuracy sought
  ///
  void Model(FitFunction *ff,
	     FitData *fd,
	     ResidualSorter *rs,
	     const int  max_iterations,
	     const REAL min_delta,
	     const REAL convergence);

};

#endif





