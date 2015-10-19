// Copyright (C) 2012-2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__BASE_GSL_HPP
#define MC__BASE_GSL_HPP

#include <iostream>
#include <sys/time.h>

#include "base_de.hpp"

#include "gsl/gsl_errno.h"
#include "gsl/gsl_odeiv2.h"

namespace mc
{
//! @brief C++ base class for computing solutions of parametric DAEs using GSL
////////////////////////////////////////////////////////////////////////
//! mc::BASE_GSL is a C++ base class for computing solutions of
//! parametric differential-algebraic equations (DAEs) using GSL
////////////////////////////////////////////////////////////////////////
class BASE_GSL: public virtual BASE_DE
{
public:
  /** @ingroup ODESLV_GSL
   *  @ingroup ODEBND_GSL
   *  @{
   */
  //! @brief Default class constructor
  BASE_GSL()
    : BASE_DE()
    {}

  //! @brief Class destructor
  virtual ~BASE_GSL()
    {}

  //! @brief GSL options
  struct Options
  {
    //! @brief Constructor
    Options():
      INTMETH(RKF45), H0(1e-2), HMIN(0e0), HMAX(0e0), NMAX(0), RTOL(1e-7), ATOL(1e-7)
      {}
    //! @brief Assignment operator
    template <typename U> Options& operator=
      ( U&options ){
        INTMETH   = options.INTMETH;
        H0        = options.H0;
        HMIN      = options.HMIN;
        HMAX      = options.HMAX;
        NMAX      = options.NMAX;
        RTOL      = options.RTOL;
        ATOL      = options.ATOL;
        return *this;
      }
    //! @brief Enumeration of numerical integration algorithms
    enum INTEGRATION_METHOD{
      RKF45=0,		//!< Explicit embedded Runge-Kutta-Fehlberg (4,5) method (non-stiff systems) [Default]
      RK8PD,		//!< Explicit embedded Runge-Kutta Prince-Dormand (8,9) method (non-stiff systems)
      MSADAMS,		//!< Variable-coefficient linear multistep Adams method in Nordsieck form (non-stiff systems)
      MSBDF		//!< Variable-coefficient linear multistep backward differentiation formula (BDF) method in Nordsieck form (stiff systems)
    };
    //! @brief Numerical integration method
    INTEGRATION_METHOD INTMETH;
    //! @brief Initial step-size (Default: 1e-2)
    double H0;
    //! @brief Minimum step-size (Default: 0e0)
    double HMIN;
    //! @brief Maximum step-size (Default: 0e0)
    double HMAX;
    //! @brief Maximum number of steps in a time stage (Default: 0)
    unsigned int NMAX;
    //! @brief Relative integration tolerance (Default: 1e-6)
    double RTOL;
    //! @brief Absolute integration tolerance (Default: 1e-6)
    double ATOL;
  };

  //! @brief Structure storing integration statistics
  struct Stats
  {
    //! @brief Constructor
    Stats():
      cputime(0.), numSteps(0), numRHS(0), numJAC(0)
      {}
    //! @brief Constructor
    void reset()
      { cputime = 0.; numSteps = numRHS = numJAC = 0; }

    //! @brief CPU time
    double cputime;
    //! @brief Number of integration steps
    unsigned int numSteps;
    //! @brief Number of right-hand side (RHS) evaluations
    unsigned int numRHS;
    //! @brief Number of Jacobian evaluations
    unsigned int numJAC;
  };
  /** @} */

protected:
  //! @brief Function to initialize GSL statistics
  static void _init_stats
    ( Stats&stats );

  //! @brief Function to finalize GSL statistics
  static void _final_stats
    ( Stats&stats );

  //! @brief Function to display GSL statistics
  static void _print_stats
    ( const Stats&stats, std::ostream&os=std::cout );

  //! @brief Private methods to block default compiler methods
  BASE_GSL(const BASE_GSL&);
  BASE_GSL& operator=(const BASE_GSL&);
};

inline void
BASE_GSL::_init_stats
( Stats&stats )
{
  // Initialize statistics
  stats.reset();
  timeval time;
  gettimeofday(&time, 0) ;
  stats.cputime = - time.tv_sec - time.tv_usec*1e-6;
}

inline void
BASE_GSL::_final_stats
( Stats&stats )
{
  // Get final CPU time
  timeval time;
  gettimeofday(&time, 0);
  stats.cputime += time.tv_sec + time.tv_usec*1e-6;
}

inline void
BASE_GSL::_print_stats
( const Stats&stats, std::ostream&os )
{
  // Statistics
  os << " No STEPS  " << stats.numSteps
     << std::endl
     << " No EVALATIONS" << "   RHS: " << stats.numRHS
                         << "   JAC: " << stats.numJAC
     << std::endl
     << " CPU TIME (SEC)     " << std::fixed << std::left
                               << std::setprecision(5) << stats.cputime
     << std::endl << std::endl;
  return;
}

} // end namescape mc

#endif

