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

#ifndef LATERATION_SORTER_HH
#define LATERATION_SORTER_HH

#include <NLMAP/IterativeModeler.hh>
#include <NLMAP/LaterationFunction.hh>

///
/// Specific multilateration-reducing
/// residual sorter for multilateration
///
class LaterationSorter : public ResidualSorter {
public:
  ///
  /// Constructor
  /// @param ff Pointer to the lateration function
  /// @param fd Pointer to the fit data
  ///
  LaterationSorter(LaterationFunction *ff, FitData *fd) : ResidualSorter(ff,fd) {}
  virtual ~LaterationSorter() {}

  ///
  /// Return maximum reesidual index
  ///
  virtual int GetMaxResidualIndex();
};

#endif
