#ifndef LEASTSQUARES_GUARD
#define LEASTSQUARES_GUARD

#include <NLMAP/Parameters.hh>
#include <NLMAP/MultiLateration.hh>

class LeastSquares {
public:
  LeastSquares(const REAL* x,
	       const REAL* y,
	       const REAL* z,
	       const REAL* d,
	       const REAL* sigma,
	       const int n);
  ~LeastSquares();

  XYZData GetPosition(const int max_it,
		      const REAL convergence);
private:
/**
   * An array of size m_number to hold the x values of each reference
   * point
   */
  REAL* m_x;

  /**
   * An array of size m_number to hold the y values of each reference
   * point
   */
  REAL* m_y;

  /**
   * An array of size m_number to hold the z values of each reference
   * point
   */
  REAL* m_z;

  /**
   * An array of size m_number to hold the distances of each reference
   * point from the located object
   */
  REAL* m_d;

  /**
   * An array of size m_number to hold the expected error for each distance
   */
  REAL* m_sigma;

  /**
   * The number of sightings for this located object
   */
  REAL m_number;
};


#endif
