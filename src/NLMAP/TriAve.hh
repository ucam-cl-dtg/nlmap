/**
 * $Header$
 */

#ifndef TRIAVE_GUARD
#define TRIAVE_GUARD

#include <NLMAP/Parameters.hh>
#include <NLMAP/MultiLateration.hh>

/**
 * An implementation of trilateration-and-average location algorithm.
 * Every combination of three readings is trilaterated and the average
 * position computed.
 */
class TriAve {
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


public:
  /**
   * Create an instance of the TriAve class which is ready to work out
   * the position for the provided data.  
   */
  TriAve(const REAL* x,
	 const REAL* y,
	 const REAL* z,
	 const REAL* d,
	 const REAL* sigma,
	 const int n);
  
  ~TriAve();
  
  /**
   * Run the Trilaterate-average algorithm and return the result if
   * successful.  Throws various exceptions (subclasses of
   * std::exception) upon failure.
   */ 
  XYZData GetPosition(const int max_it,
		      const REAL convergence);

private:
  
  /**
   * Choose all the permutations by fixing the first choice and then
   * calling Choose2 for each one.
   */
  void TriAve::Choose1(REAL& sumx, 
		       REAL& sumy, 
		       REAL& sumz, 
		       REAL& sumxsq,
		       REAL& sumysq,
		       REAL& sumzsq,
		       int& count,
		       const int max_it,
		       const REAL convergence);

  /**
   * Choose all the permutations by choosing the second choice
   * (excluding that chosen for the first choice) and then calling
   * Choose3 for each one.
   */
  void TriAve::Choose2(int first,
		       REAL& sumx, 
		       REAL& sumy, 
		       REAL& sumz, 
		       REAL& sumxsq,
		       REAL& sumysq,
		       REAL& sumzsq,
		       int& count,
		       const int max_it,
		       const REAL convergence);
  
  /**
   * Choose all the permutations by choosing the third choice
   * (excluding the first and second selections already made) and then
   * calling Process for each one.
   */
  void TriAve::Choose3(int first, 
		       int second,
		       REAL& sumx, 
		       REAL& sumy, 
		       REAL& sumz, 
		       REAL& sumxsq,
		       REAL& sumysq,
		       REAL& sumzsq,
		       int& count,
		       const int max_it,
		       const REAL convergence);


  /**
   * Trilaterate the chosen triplet and accumulate the averages.
   */
  void TriAve::Process(int first, 
		       int second, 
		       int third,
		       REAL& sumx, 
		       REAL& sumy, 
		       REAL& sumz, 
		       REAL& sumxsq,
		       REAL& sumysq,
		       REAL& sumzsq,
		       int& count,
		       const int max_it,
		       const REAL convergence);
};

#endif//TRIAVE_GUARD
