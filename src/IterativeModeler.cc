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

#include <NLMAP/IterativeModeler.hh>
#include <NLMAP/NonLinearModel.hh>
#include <iostream>



ResidualSorter::ResidualSorter(FitFunction *ff, FitData *fd)
  : mFunc(ff), mData(fd) {}

ResidualSorter::~ResidualSorter() {};

void IterativeModeler::Model(FitFunction *ff,
			     FitData *fd,
			     ResidualSorter *rs,
			     const int  max_it,
			     const REAL min_delta,
			     const REAL convergence)
{
  ff->InitialiseParameters(fd);

  if (fd->GetInputDataSize() <= ff->GetNumParameters()) {
    // Not enough constraints
    throw ModelingError("Too few constraints");
  }
  
  bool success=false;
  while (fd->GetInputDataSize()>ff->GetNumParameters()) {
    // Try to fit a model to the current parameters
    ff->InitialiseParameters(fd);
    NonLinearModel nlm(ff,fd);
    try {
      nlm.Fit(max_it, min_delta); 
    }
    catch (NLMAPException &nlme) {
      throw ModelingError("NLM exception");
    }
    ff->SetError(nlm.GetStdErr());
    // Is the result accurate enough?
    if (nlm.GetStdErr() < convergence) {
      success=true;
      break;
    }
    else {
      // Need to discard data and try again
      try {
	int n = rs->GetMaxResidualIndex();
	fd->Invalidate(n);
      }
      catch(...) { throw ModelingError("Failed to ID max residual"); }
    }
  }

  if (!success) throw ModelingFailure();
}
