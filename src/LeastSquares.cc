
#include <NLMAP/LeastSquares.hh>
#include <NLMAP/NonLinearModel.hh>
#include <NLMAP/LaterationFunction.hh>
#include <cmath>
#include <cstdlib>

LeastSquares::LeastSquares(const REAL* x,
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

LeastSquares::~LeastSquares() {
  delete[] m_x;
  delete[] m_y;
  delete[] m_z;
  delete[] m_d;
  delete[] m_sigma;
}


XYZData LeastSquares::GetPosition(const int max_it,
				  const REAL convergence)
{

  LaterationData latdata;
  for(int i=0;i<m_number;++i) {
    LaterationData::Datum d;
    d.push_back(m_x[i]);
    d.push_back(m_y[i]);
    d.push_back(m_z[i]);
    latdata.AddDatum(d,m_d[i],m_sigma[i]);
  }
  LaterationFunction ff;
  try {
    ff.InitialiseParameters(&latdata);
    NonLinearModel nlm(&ff,&latdata);
    nlm.Fit(max_it, convergence);

    XYZData xyzd;
    xyzd.x = ff.GetParams()[0];
    xyzd.y = ff.GetParams()[1];
    xyzd.z = ff.GetParams()[2];
    xyzd.sigma = nlm.GetStdErr();
    return xyzd;
  }
  catch (ModelingFailure& e) {
    DEBUG("Caught ModelingFailure - skipping this one");
  }
  catch (MaxIterations& e) {
    DEBUG("Caught MaxIterations - skipping this one");
  }  
}
