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

#ifndef ITERATIVE_MODELER_HH
#define ITERATIVE_MODELER_HH

#include <NLMAP/Parameters.hh>
#include <NLMAP/NonLinearModel.hh>
#include <NLMAP/NLMAPExceptions.hh>
#include <cstring>

/**
 * Pure virtual class to override with your
 * specifics for ordering the model residuals
 * and discarding outliers
 */
class ResidualSorter {
public:
  ResidualSorter(FitFunction *ff, FitData *fd);
  virtual ~ResidualSorter();
  // Return the index to discard
  virtual int GetMaxResidualIndex()=0;
protected:
  FitFunction  *mFunc;
  FitData      *mData;
};




class IterativeModeler {
public:
  IterativeModeler() {};

  void Model(FitFunction *ff,
	     FitData *fd,
	     ResidualSorter *rs,
	     const int  max_iterations,
	     const REAL min_delta,
	     const REAL convergence);

protected:
};

#endif





