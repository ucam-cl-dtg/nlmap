/**
 * $Header$
 */

#include <NLMAP/TriAve.hh>

#include <iostream>
#include <cmath>
#include <cstdlib>

#define MAX_ERROR 0.05

void do_attempt(REAL tx, REAL ty, REAL tz,TriAve& r1) {
  try {
    XYZData pd = r1.GetPosition(100,0.001);
    REAL distance = sqrt((pd.x - tx)*(pd.x - tx) + (pd.y-ty)*(pd.y-ty) + (pd.z-tz)*(pd.z-tz));
    
    std::cout << "Distance was " << distance << std::endl;
  }
  catch (ModelingFailure &mf) {
    std::cout << "Caught modelling failure - " << mf.what() << " - test failed."<< std::endl;        
    exit(-1);
  }
  catch (NLMAPException &n) {
    std::cerr << n.what() << std::endl;
    exit(-1);
  }

  std::cout << "Pass" << std::endl;
}

REAL peturb(REAL value) {
  return value + (value*0.00001*(rand() - RAND_MAX/2))/(RAND_MAX/2);  
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
  

  //  for(int i=0;i<9;i++) {
  //    std::cout << "Point " << i << " ("<<x[i]<<","<<y[i]<<","<<z[i]<<","<<d[i]<<")"<<std::endl;
  //  }

  // try it first with three points
  for(int i=0;i<7;++i) {
    std::cout << "Running with points " << i << "," << i+1 << "," << i+2 << std::endl;
    TriAve r1(x+i,y+i,z+i,d+i,s+i,3);
    do_attempt(tx,ty,tz,r1);
  }
  

  for(int i=0;i<9;++i) {
    // add an outrageous outlier
    REAL store = d[i];
    d[i] *= 1.5;
    std::cout << "Running with all points and one artificial outlier at point " << i << std::endl;
    TriAve r1(x,y,z,d,s,9);
    do_attempt(tx,ty,tz,r1);
    d[i] = store;
  }

  for(int i=0;i<9;++i) {
    // add an outrageous outlier
    REAL store1 = d[i];
    d[i] *= 1.5;
    for(int j=0;j<9;++j) {
      if (i != j) {
	REAL store2 = d[j];
	d[j] *= 1.5;
	std::cout << "Running with all points and artificial outliers at points " << i << " and " << j << std::endl;
	TriAve r1(x,y,z,d,s,9);
	do_attempt(tx,ty,tz,r1);
	d[j] = store2;
      }
      d[i] = store1;
    }
  }
}

/**
 * Test harness for the TriAve implementation
 */
int main(int argc, char **argv) {

  test1(1,1,1,1);
  //  test1(1,1,1,2);
  //  test1(-1,5,-45,13);

}
