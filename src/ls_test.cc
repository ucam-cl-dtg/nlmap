#include <NLMAP/LinearisedLS.hh>
#include <cmath>
#include <NLMAP/MultiLateration.hh>

#include <iostream>
int main(int argc, char **argv) {
 // Set up for somewhere near (1,1,0)
  REAL x[5] = {0.0, 2.0, 0.0, 1.0, 0.0};
  REAL y[5] = {0.0, 0.0, 2.0, 1.0, 0.0};
  REAL z[5] = {1.0, 1.0, 1.0, 2.0, 1.0};
  // REAL d[4] = {sqrt(3.0),sqrt(3.0),sqrt(3.0),2.0};
  REAL d[5] = {1.7, 1.73, 1.71, 2.1, 21.3};
  REAL s[5] = {1.0, 1.0, 1.0, 1.0, 10.0};

  LinearisedLS ls(x,y,z,d,s,4);

  XYZData pd = ls.GetPosition();
  
  std::cout << pd.x << " " << pd.y << " " << pd.z << std::endl;

  MultiLateration ml(x,y,z,d,s,5);
 try {
  XYZData pd2 = ml.GetPosition(100,0.001,0.03);
 }
 catch (ModelingFailure &mf) {
   std::cerr << mf.what() << std::endl;
   std::cout << "Best attempt: " << "( " << ml.GetLaterationFunction()->GetParams()[0]
	      << ", " <<ml.GetLaterationFunction()->GetParams()[1] << "," <<
      ml.GetLaterationFunction()->GetParams()[2] << ")" << std::endl;
    std::cout << "With error: " << ml.GetLaterationFunction()->GetError() << std::endl;
 }


  
};
