#include <NLMAP/Parameters.hh>
#include <NLMAP/MultiLateration.hh>

class RobustLS {

public:
  RobustLS(REAL* x,
	   REAL* y,
	   REAL* z,
	   REAL* d,
	   REAL* sigma,
	   int n);
  ~RobustLS();

  XYZData GetPosition(const int maxit);

private:
  void LeastSquares();

  REAL **mA;
  REAL  *mB;
  REAL **mH;
  REAL  *mS;
  REAL   mP[3];
  int    mDim;
  int    mN;
  
  REAL  *mX;
  REAL  *mY;
  REAL  *mZ;

};
