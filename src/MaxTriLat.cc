/*
  Copyright (C) 2004 Andrew C. Rice

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

  Email: acr31@cam.ac.uk
*/

/**
 * $Header$
 */


#include <NLMAP/MaxTriLat.hh>

MaxTriLat::MaxTriLat(const REAL* x,
		     const REAL* y,
		     const REAL* z,
		     const REAL* d,
		     const REAL* sigma,
		     const int n) : MinTriLat(x,y,z,d,sigma,n)
{}

void MaxTriLat::MakeSelection(int ids[3]) const {
  for(int i=0;i<m_number;++i) {
    if (m_d[i] > m_d[ids[0]]) {
      ids[2] = ids[1];
      ids[1] = ids[0];
      ids[0] = i;
    }
    else if (m_d[i] > m_d[ids[1]]) {
      ids[2] = ids[1];
      ids[1] = i;
    }
    else if (m_d[i] > m_d[ids[2]]) {
      ids[2] = i;
    }
  }
}
