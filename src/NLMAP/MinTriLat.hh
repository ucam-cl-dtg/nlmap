/*
  Copyright (C) 2004 Andrew C. Rice

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

  Email: acr31@cam.ac.uk
*/

/**
 * $Header$
 */

#ifndef MINTRILAT_GUARD
#define MINTRILAT_GUARD

#include <NLMAP/Parameters.hh>
#include <NLMAP/MultiLateration.hh>

/**
 * An implementation of minimal tri-lateration.  The three shortest
 * readings are tri-laterated and the position returned.
 */
class MinTriLat {
protected:
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

  virtual void MakeSelection(int ids[3]) const;

public:
  /**
   * Create an instance of the MinTriLat class which is ready to work out
   * the position for the provided data.  
   */
  MinTriLat(const REAL* x,
	    const REAL* y,
	    const REAL* z,
	    const REAL* d,
	    const REAL* sigma,
	    const int n);
  
  virtual ~MinTriLat();
  
  /**
   * Run the minimal Trilaterate algorithm and return the result if
   * successful.  Throws various exceptions (subclasses of
   * std::exception) upon failure.
   */ 
  XYZData GetPosition(const int max_it,
		      const REAL convergence);
  
};

#endif//MINTRILAT_GUARD
