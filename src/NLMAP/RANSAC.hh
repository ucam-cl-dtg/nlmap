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

#ifndef RANSAC_GUARD
#define RANSAC_GUARD

#include <NLMAP/Parameters.hh>
#include <NLMAP/MultiLateration.hh>


/**
 * An implementation of the RANdom SAmple Consensus algorithm.  See:
 * "Random Sample Consensus: A Paradigm for Model Fitting with
 * Applications to Image Analysis and Automated Cartography", in
 * Communications of the ACM, June 1981, Vol 24, Number 6.  Also see:
 * http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/FISHER/RANSAC/
 */
class RANSAC {
private:
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

  /**
   * The probability of a randomly selected data item being part of a
   * good model
   */
  REAL m_pgood;

  /**
   * The probability that the algorithm will exit without finding a
   * good fit if one exists (false negative)
   */
  REAL m_pfail;

  /**
   * The number of times to search for a valid quorum of data
   * points. This is calculated from p_fail and p_good in the
   * following manner: pfail = p(L consecutive failures) pfail =
   * p(given trial is a failure)**L pfail = p(1-given trial
   * succeeds)**L pfail = p(1-p_good**N)**L L =
   * log(p_fail)/log(1-p_good**N) Note: N is the number of data points
   * required for tri-lateration (three).  L is the number of
   * iterations to perform before giving up.
   */
  int m_iteration_count;

public:
  /**
   * Create an instance of the RANSAC class that is prepared to
   * calculate a position from the provided sightings.  The data in
   * arrays x,y,z,d and sigma are copied. pgood (the probability of a
   * randomly selected data item being part of a good model) and pfail
   * (the probability that the algorithm will exit without finding a
   * good fit if one exists) are used to calculate the number of
   * retries to be made before giving up.
   */ 
  RANSAC(const REAL* x,
	 const REAL* y,
	 const REAL* z,
	 const REAL* d,
	 const REAL* sigma,
	 const int n,
	 const REAL pgood,
	 const REAL pfail);

  ~RANSAC();

  /**
   * Run the RANSAC algorithm and return the result if successful.
   * Throws various exceptions (subclasses of std::exception) upon
   * failure.
   */ 
  XYZData GetPosition(const int max_it,
		      const REAL convergence);
};

#endif//RANSAC_GUARD
