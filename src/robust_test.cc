#include <NLMAP/RobustLS.hh>
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
  REAL s[5] = {0.4,0.4,0.4,0.4,0.4};

  RobustLS rs(x,y,z,d,s,5);

  XYZData pd = rs.GetPosition(100);
  
  std::cout << pd.x << " " << pd.y << " " << pd.z << std::endl;
  
};
