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

#ifndef POSITIONER_HH
#define POSITIONER_HH

#include <NLMAP/Parameters.hh>
#include <NLMAP/LaterationFunction.hh>
#include <NLMAP/LaterationSorter.hh>
#include <NLMAP/AngulationFunction.hh>
#include <NLMAP/MagnitudeSorter.hh>
#include <exception>
#include <cstring>

struct PositionData {
  REAL x;
  REAL y;
  REAL z;
  REAL sigma;
};


class PositioningError {
public:
  PositioningError(char *f){strcpy(msg,f);}
  PositioningError() {strcpy(msg,"Error in positioning algorithm");}
  char * what() { return msg;}
private: 
  char  msg[512];
};

class Positioner {
public:
  Positioner() {};
  PositionData & GetPosition(LaterationFunction *lf,
			     LaterationSorter *ls,
			     const int  max_it,
			     const REAL min_delta,
			     const REAL convergence);

  PositionData & GetPosition(AngulationFunction *lf,
			     MagnitudeSorter *ls,
			     const int  max_it,
			     const REAL min_delta,
			     const REAL convergence);
  
protected:
  void GetPositionEstimate(LaterationFunction *lf); 
  PositionData         mResult;
};

#endif





