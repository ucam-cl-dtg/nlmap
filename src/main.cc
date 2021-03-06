/*
  Copyright (C) 2004 Robert K. Harle

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Email: robert.harle@achilles.org
  Email: rkh23@cantab.net
*/


#include <NLMAP/MultiLateration.hh>
#include <NLMAP/MultiAngulation.hh>

#include <iostream>
#include <cmath>

//----------------------------------
// This file is an example of how to use
// the non linear modeller
//----------------------------------
int main(int argc, char **argv) {

  std::cout << "BAT: 996.949 1041.04 0.704703" << std::endl;
 



	REAL x[9] = {998.379822,998.379700,997.556274,997.107727,997.105835,996.378906,995.789734,994.583435,995.186523};
  REAL y[9] = {1042.117554,1040.692261,1041.737427,1040.835205,1042.340332,1041.953247,1041.354980,1040.754395,1042.480225};
  REAL z[9] = {2.282333,2.278333,2.282333,2.278333,2.286333,2.467333,2.463667,2.470667,2.465667};
  REAL d[9] = {2.391647,2.146659,1.823898,1.600946,2.039072,2.070182,2.122031,3.521961,2.879030};
  REAL s[9] = {1.0, 1.0, 1.0, 1.0, 1.0,1.0,1.0,1.0,1.0};




 MultiLateration ml(x,y,z,d,s,9);
 try {
   XYZData pd = ml.GetPosition(100,0.001,0.03);
 }
 catch (ModelingFailure &mf) {
   std::cerr << mf.what() << std::endl;
   std::cout << "Best attempt: " << "( " << ml.GetLaterationFunction()->GetParams()[0]
	      << ", " <<ml.GetLaterationFunction()->GetParams()[1] << "," <<
      ml.GetLaterationFunction()->GetParams()[2] << ")" << std::endl;
    std::cout << "With error: " << ml.GetLaterationFunction()->GetError() << std::endl;
 }
 return 0;



  //----------------------------------
  // MULTILATERATION
  // First we test multilateration, where
  // we have a series of distance ranges to bases
  // at known locations and want to know
  // the most likely source position.
  //----------------------------------

  // Set up for somewhere near (1,1,0)
 //   REAL x[5] = {0.0, 2.0, 0.0, 1.0, 0.0};
//    REAL y[5] = {0.0, 0.0, 2.0, 1.0, 0.0};
//    REAL z[5] = {1.0, 1.0, 1.0, 2.0, 1.0};
//    REAL d[5] = {1.7, sqrt(3.0), sqrt(3.0), 2.1, sqrt(3.0)};
//    REAL s[5] = {1.0, 1.0, 1.0, 1.0, 10.0};

 // MultiLateration ml(x,y,z,d,s,5);
 
  try {
    XYZData pd = ml.GetPosition(100,0.001,0.03);
    std::cout << "*** Multilateration" << std::endl;
    std::cout << "Result (should be close to (1,1,0): " <<
      "(" << pd.x << " " << pd.y << " " << pd.z << ") " << std::endl;
    std::cout << "Standard error: " << pd.sigma << " distance units" << std::endl;
  
  }
  catch (ModelingFailure &mf) {
    std::cerr << mf.what() << std::endl;
    std::cout << "Best attempt: " << "( " << ml.GetLaterationFunction()->GetParams()[0]
	 << ", " <<ml.GetLaterationFunction()->GetParams()[1] << "," <<
      ml.GetLaterationFunction()->GetParams()[2] << ")" << std::endl;
    std::cout << "With error: " << ml.GetLaterationFunction()->GetError() << std::endl;
  }
  catch (NLMAPException &n) {
    std::cerr << n.what() << std::endl;
  }

  std::cout << std::endl;

  //----------------------------------
  // MULTIANGULATION
  // Now we test multiangulation, where
  // we have a series of bases
  // at known locations, each with a cw bearing
  // from (0,1) to a source at unknown (x,y)
  //----------------------------------
  
  REAL x2[5] = {0.0, 1.0, 1.0, 1.0};
  REAL y2[5] = {2.0, 1.0, 0.0, 2.0};
  REAL b[5] = {181.0/180*M_PI,
	       271.0/180*M_PI,
	       316.0/180*M_PI,
	       302.0/180*M_PI};
  REAL s2[5] = {10.0/180*M_PI,
		10.0/180*M_PI,
		10.0/180*M_PI,
		10.0/180*M_PI};

  std::cout << "*** Multiangulation" << std::endl;

  XYData xyd;
  MultiAngulation ma(x2,y2,b,s2,4);
  try {
    xyd = ma.GetPosition(100,0.001,40/180*M_PI);
    
    std::cout << "Result (should be close to (0,1): " <<
      "(" << xyd.x << " " << xyd.y << ") " << std::endl;
    std::cout << "Standard error: " << xyd.sigma << " radians" << std::endl;
  }
  catch (ModelingFailure &mf) {
    std::cerr << mf.what() << std::endl;
    std::cout << "Best attempt: " << "( " << ma.GetAngulationFunction()->GetParams()[0]
	 << ", " <<ma.GetAngulationFunction()->GetParams()[1] << ")" << std::endl;
    std::cout << "With error: " << ma.GetAngulationFunction()->GetError() << std::endl;
  }
  catch (NLMAPException &n) {
    std::cerr << n.what() << std::endl;
  }

}
