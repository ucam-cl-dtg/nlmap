/**
 * $Header$
 */

#include <NLMAP/RANSAC.hh>
#include <NLMAP/NonLinearModel.hh>
#include <NLMAP/LaterationFunction.hh>
#include <cmath>
#include <cstdlib>

RANSAC::RANSAC(const REAL* x,
	       const REAL* y,
	       const REAL* z,
	       const REAL* d,
	       const REAL* sigma,
	       const int n,
	       const REAL pgood,
	       const REAL pfail): 
  m_x(new REAL[n]),
  m_y(new REAL[n]),
  m_z(new REAL[n]),
  m_d(new REAL[n]),
  m_sigma(new REAL[n]),
  m_number(n),
  m_pgood(pgood),
  m_pfail(pfail)  
{
  DEBUG("Pgood: " << pgood);
  DEBUG("Pfail: " << pfail);
  
  for(int i=0;i<n;++i) {
    m_x[i] = x[i];
    m_y[i] = y[i];
    m_z[i] = z[i];
    m_d[i] = d[i];
    m_sigma[i] = sigma[i];
  }

  m_iteration_count = (int)(log(pfail)/log(1-pgood*pgood*pgood));

  DEBUG("Iteration count: " << m_iteration_count);
}

RANSAC::~RANSAC() {
  delete[] m_x;
  delete[] m_y;
  delete[] m_z;
  delete[] m_d;
  delete[] m_sigma;
}

XYZData RANSAC::GetPosition(const int max_it,
			    const REAL convergence)
{
  
  for(int iterations=0;iterations<m_iteration_count;++iterations) {

    // select 3 data points at random making sure we have no
    // collisions
    int points[3];
    for(int i=0;i<3;++i) {
      points[i] = (int)((REAL)(m_number*rand())/RAND_MAX);
      for(int j=0;j<i;) {
	if (points[i] == points[j]) {
	  points[i] = (int)((REAL)(m_number*rand())/RAND_MAX);	  
	  j=0;
	}
	else {
	  ++j;
	}
      }
    }

    DEBUG("Points selection: [" << points[0] << "," << points[1] << "," << points[2] << "]");

    // load them into the Nonlinear lateration fit function
    LaterationData latdata;
    for(int i=0;i<3;++i) {
      int index = points[i];
      LaterationData::Datum d;
      d.push_back(m_x[index]);
      d.push_back(m_y[index]);
      d.push_back(m_z[index]);
      latdata.AddDatum(d,m_d[index],m_sigma[index]);
    }
    LaterationFunction ff(3,&latdata);
    ff.InitialiseParameters();
    NonLinearModel nlm(&ff,&latdata);
    nlm.Fit(max_it, convergence); 
         
    REAL nlm_x = ff.GetParams()[0];
    REAL nlm_y = ff.GetParams()[1];
    REAL nlm_z = ff.GetParams()[2];

    // this is only one of two possible positions the other one lies
    // on the other sied of the plane described by these three points.

    // calculate the normal vector to the plane by working out two
    // vectors that lie on the plane and taking their cross product

    REAL v1x = m_x[points[0]] - m_x[points[1]];
    REAL v1y = m_y[points[0]] - m_y[points[1]];
    REAL v1z = m_z[points[0]] - m_z[points[1]];

    DEBUG("v1 = (" << v1x << "," << v1y << "," << v1z << ")");

    REAL v2x = m_x[points[0]] - m_x[points[2]];
    REAL v2y = m_y[points[0]] - m_y[points[2]];
    REAL v2z = m_z[points[0]] - m_z[points[2]];

    DEBUG("v2 = (" << v2x << "," << v2y << "," << v2z << ")");

    REAL nx = v1y*v2z - v2y*v1z;
    REAL ny = v2x*v1z - v1x*v2z;
    REAL nz = v1x*v2y - v1y*v2x;

    REAL len = sqrt( nx*nx + ny*ny + nz*nz);
    
    nx/=len;
    ny/=len;
    nz/=len;

    DEBUG("Normal vector to plane is " << nx << "," << ny << "," << nz);

    // project the vector from the located point to the plane onto the
    // normal vector
    REAL dotprod = nx*(nlm_x - m_x[points[0]]) +  ny*(nlm_y - m_y[points[0]]) +  nz*(nlm_z - m_z[points[0]]);
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


    // find out how many points within the data set fit the
    // model
    int t =0;
    REAL sum = 0;
    REAL sumsq = 0;
    for(int i=0;i<m_number;++i) {
      REAL dist = sqrt( (m_x[i]-nlm_x)*(m_x[i]-nlm_x) + 
			(m_y[i]-nlm_y)*(m_y[i]-nlm_y) + 
			(m_z[i]-nlm_z)*(m_z[i]-nlm_z) );
      REAL error = fabs( m_d[i] - dist );
      
      DEBUG("Sighting " << i << " error: " << error);
      sum += error;
      sumsq += error*error;

      if (error < m_sigma[i]) {
	t++;
      }
    }

    DEBUG("t: " << t);

    // if this is more than 8 points (see the original paper for where
    // this number comes from) then return the data
    if (t == m_number || t > 8) {
      XYZData pd;
      pd.x = nlm_x;
      pd.y = nlm_y;
      pd.z = nlm_z;
      pd.sigma = sqrt((sum*sum - sumsq)/t);
      return pd;
    }

    // if this is less than 8 points then repeat
    
  }
  
  // if we get here we don't deem it likely that we will get a good
  // answer - so give up.
  throw ModelingFailure();
  

}
