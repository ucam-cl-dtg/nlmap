
#include <NLMAP/Parameters.hh>
#include <NLMAP/MultiLateration.hh>

///
/// Perform multilateration by using 
/// Linearised Least Squares
///
class LinearisedLS {

public:

  LinearisedLS(const REAL* x,
	       const REAL* y,
	       const REAL* z,
	       const REAL* d,
	       const REAL* sigma,
	       const int n);

  virtual ~LinearisedLS();

  ///
  /// Get the position using linearised
  /// least squares
  ///
  XYZData GetPosition();

  ///
  /// Get the hat matrix diagonals
  ///
  REAL * GetHatMatrixDiagonals() { return mH;}

private:
  REAL **mA;
  REAL  *mB;
  REAL  *mH;
  int mDim;

};
