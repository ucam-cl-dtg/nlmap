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

#include <NLMAP/MinTriLat.hh>
#include <cmath>

MinTriLat::MinTriLat(const REAL* x,
	       const REAL* y,
	       const REAL* z,
	       const REAL* d,
	       const REAL* sigma,
	       const int n) :
  m_x(new REAL[n]),
  m_y(new REAL[n]),
  m_z(new REAL[n]),
  m_d(new REAL[n]),
  m_sigma(new REAL[n]),
  m_number(n)
{

  for(int i=0;i<n;++i) {
    m_x[i] = x[i];
    m_y[i] = y[i];
    m_z[i] = z[i];
    m_d[i] = d[i];
    m_sigma[i] = sigma[i];
  }
}

MinTriLat::~MinTriLat() {
  delete[] m_x;
  delete[] m_y;
  delete[] m_z;
  delete[] m_d;
  delete[] m_sigma;
}

void MinTriLat::MakeSelection(int ids[3]) const {
  for(int i=0;i<m_number;++i) {
    if (m_d[i] < m_d[ids[0]]) {
      ids[2] = ids[1];
      ids[1] = ids[0];
      ids[0] = i;
    }
    else if (m_d[i] < m_d[ids[1]]) {
      ids[2] = ids[1];
      ids[1] = i;
    }
    else if (m_d[i] < m_d[ids[2]]) {
      ids[2] = i;
    }
  }
}

XYZData MinTriLat::GetPosition(const int max_it,
			       const REAL convergence) {
  // take the three sightings correlating to the smallest distances
  int smallest[3] = {0};
  MakeSelection(smallest);
  
  const int first = smallest[0];
  const int second  = smallest[1];
  const int third = smallest[2];
  
  DEBUG("Processing " << first << "," << second << "," << third);

  // calculate the normal vector to the plane by working out two
  // vectors that lie on the plane and taking their cross product

  REAL v1x = m_x[first] - m_x[second];
  REAL v1y = m_y[first] - m_y[second];
  REAL v1z = m_z[first] - m_z[second];
  
  DEBUG("v1 = (" << v1x << "," << v1y << "," << v1z << ")");
  
  REAL v2x = m_x[first] - m_x[third];
  REAL v2y = m_y[first] - m_y[third];
  REAL v2z = m_z[first] - m_z[third];
  
  DEBUG("v2 = (" << v2x << "," << v2y << "," << v2z << ")");
  
  REAL nx = v1y*v2z - v2y*v1z;
  REAL ny = v2x*v1z - v1x*v2z;
  REAL nz = v1x*v2y - v1y*v2x;
  
  DEBUG("n = (" << nx << "," << ny << "," << nz << ")" );
  
  // compute the magnitude of the cross product - if this is close
  // to zero then our vectors are close to collinear and so we
  // should not consider them
  
  REAL mag = nx*nx+ny*ny+nz*nz;
  
  if (mag < 0.0001) {
    DEBUG("Magnitude of cross product is too small - the selected vectors are collinear - bailing out.");
    throw ModelingError("Magnitude of cross product is too small");
  }
  
  REAL len = sqrt(mag);
  
  nx/=len;
  ny/=len;
  nz/=len;
  
  DEBUG("Normal vector to plane is " << nx << "," << ny << "," << nz);
  
  
  // load the chosen points into the NonLinear lateration function  
  LaterationData latdata;

  LaterationData::Datum d1;
  d1.push_back(m_x[first]);
  d1.push_back(m_y[first]);
  d1.push_back(m_z[first]);
  latdata.AddDatum(d1,m_d[first],m_sigma[first]);

  LaterationData::Datum d2;
  d2.push_back(m_x[second]);
  d2.push_back(m_y[second]);
  d2.push_back(m_z[second]);
  latdata.AddDatum(d2,m_d[second],m_sigma[second]);

  LaterationData::Datum d3;
  d3.push_back(m_x[third]);
  d3.push_back(m_y[third]);
  d3.push_back(m_z[third]);
  latdata.AddDatum(d3,m_d[third],m_sigma[third]);

  LaterationFunction ff;
  ff.InitialiseParameters(&latdata);
  NonLinearModel nlm(&ff,&latdata);
  nlm.Fit(max_it, convergence); 

  REAL nlm_x = ff.GetParams()[0];
  REAL nlm_y = ff.GetParams()[1];
  REAL nlm_z = ff.GetParams()[2];
  
  // project the vector from the located point to the plane onto the
  // normal vector
  REAL dotprod = nx*(nlm_x - m_x[first]) +  ny*(nlm_y - m_y[first]) +  nz*(nlm_z - m_z[first]);
  REAL px = nx*dotprod;
  REAL py = ny*dotprod;
  REAL pz = nz*dotprod;
  
  // now work out the new point by adding twice this vector onto the located point
  REAL nlm_x2 = nlm_x-2*px;
  REAL nlm_y2 = nlm_y-2*py;
  REAL nlm_z2 = nlm_z-2*pz;
  
  DEBUG("Resulting position 1 is ("<<nlm_x<<","<<nlm_y<<","<<nlm_z<<")");
  DEBUG("Resulting position 2 is ("<<nlm_x2<<","<<nlm_y2<<","<<nlm_z2<<")");
  
  if (nlm_z2 < nlm_z) {
    DEBUG("Swapping position (selecting minumum z value");
    nlm_x = nlm_x2;
    nlm_y = nlm_y2;
    nlm_z = nlm_z2;
  }
  
  XYZData result;
  result.x = nlm_x;
  result.y = nlm_y;
  result.z = nlm_z;
  result.sigma = 100;
  return result;
}

