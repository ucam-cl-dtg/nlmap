/**
 * $Header$
 */

#include <NLMAP/RANSAC.hh>

#include <iostream>
#include <cmath>
#include <cstdlib>

#define MAX_ERROR 0.01

void do_attempt(REAL tx, REAL ty, REAL tz,RANSAC& r1) {
  try {
    XYZData pd = r1.GetPosition(100,0.001);
    REAL distance = sqrt((pd.x - tx)*(pd.x - tx) + (pd.y-ty)*(pd.y-ty) + (pd.z-tz)*(pd.z-tz));
    
    if (distance > MAX_ERROR) {
      std::cout << "Distance was " << distance << " which exceeds threshold ("<<MAX_ERROR<<") - test failed." << std::endl;
      exit(-1);
    }  
  }
  catch (ModelingFailure &mf) {
    std::cout << "Caught modelling failure - " << mf.what() << " - test failed."<< std::endl;


    exit(-1);
  }
  catch (NLMAPException &n) {
    std::cerr << n.what() << std::endl;
    exit(-1);
  }
}

REAL peturb(REAL value) {
  return value + (value*0.001*(rand() - RAND_MAX/2))/(RAND_MAX/2);  
}

void test1(REAL tx, REAL ty, REAL tz, REAL d1) {
  std::cout << "Running test1 around target at ("<<tx<<","<<ty<<","<<tz<<")"<<std::endl;

  // these points are a 3x3 grid centred directly above the target point
  // they are cunningly ordered so that any group of three are non-collinear
  
  REAL x[] = {tx-d1,tx-d1,tx,   tx+d1,tx-d1,tx,   tx+d1,tx,   tx+d1};
  REAL y[] = {ty-d1,ty+d1,ty-d1,ty-d1,ty,   ty+d1,ty,   ty,   ty+d1};
  REAL z[] = {tz+d1,tz+d1,tz+d1,tz+d1,tz+d1,tz+d1,tz+d1,tz+d1,tz+d1};
  REAL d[9];
  REAL s[9];
  
  for(int i=0;i<9;++i) {
    d[i] = peturb(sqrt( (x[i]-tx)*(x[i]-tx)+(y[i]-ty)*(y[i]-ty)+(z[i]-tz)*(z[i]-tz) ));
    s[i] = 0.1;
  }
  

  for(int i=0;i<9;i++) {
    std::cout << "Point " << i << " ("<<x[i]<<","<<y[i]<<","<<z[i]<<","<<d[i]<<")"<<std::endl;
  }

  // try it first with three points
  for(int i=0;i<7;++i) {
    std::cout << "Running with points " << i << "," << i+1 << "," << i+2 << std::endl;
    RANSAC r1(x+i,y+i,z+i,d+i,s+i,3,0.99,0.001);
    do_attempt(tx,ty,tz,r1);
  }

}

/**
 * Test harness for the RANSAC implementation
 */
int main(int argc, char **argv) {

  test1(1,1,1,1);

}
