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

#include <NLMAP/Positioner.hh>

#include <iostream>
#include <cmath>

//----------------------------------
// This file is an example of how to use
// the non linear modeller
//----------------------------------
int main(int argc, char **argv) {

  LaterationData data;
  FitData::Datum datum;
  datum.push_back(0.0);
  datum.push_back(0.0);
  datum.push_back(2.0);
  data.AddDatum(datum,1.8499,0.1);

  datum.clear();
  datum.push_back(1.0);
  datum.push_back(0.0);
  datum.push_back(2.0);
  data.AddDatum(datum,1.38434,0.1);

  datum.clear();
  datum.push_back(2.0);
  datum.push_back(0.0);
  datum.push_back(2.0);
  data.AddDatum(datum,1.83012,0.1);


  LaterationFunction lf(3,&data);
  lf.InitialiseParameters();
  NonLinearModel nlm(&lf,&data);
  try{
    nlm.Fit(100, 0.001); 
    std::cout << "Position: " << lf.GetParams()[0] << " " << lf.GetParams()[1] << " " << lf.GetParams()[2] << std::endl;
  }
  catch (NLMAPException &n) {
    std::cout << n.what() << std::endl;
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
  REAL x[5] = {0.0, 2.0, 0.0, 1.0, 0.0};
  REAL y[5] = {0.0, 0.0, 2.0, 1.0, 0.0};
  REAL z[5] = {1.0, 1.0, 1.0, 2.0, 1.0};
  REAL d[5] = {1.7, sqrt(3.0), sqrt(3.0), 2.1, sqrt(3.0)};
  REAL s[5] = {1.0, 1.0, 1.0, 1.0, 10.0};

  MultiLateration ml(x,y,z,d,s,5);
 
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
