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
};


#endif
