// Copyright (C) 2015 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__BASE_EXPAND_HPP
#define MC__BASE_EXPAND_HPP

#include <iostream>
#include <iomanip>
#include <sys/time.h>

#include "mcfunc.hpp"
#include "base_rk.hpp"

namespace mc
{
//! @brief C++ base class for computing solutions of parametric DAEs using GSL
////////////////////////////////////////////////////////////////////////
//! mc::BASE_EXPAND is a C++ base class for computing solutions of
//! parametric differential-algebraic equations (DAEs) using
//! discretization/expansion techniques
////////////////////////////////////////////////////////////////////////
class BASE_EXPAND:
  protected BASE_RK
{
public:
  //! @brief Default class constructor
  BASE_EXPAND()
    {}

  //! @brief Class destructor
  virtual ~BASE_EXPAND()
    {}

  //! @brief GSL options
  struct Options
  {
    //! @brief Constructor
    Options():
      INTMETH(TS), H0(1e-2), HMIN(0e0), HMAX(0e0), RTOL(1e-8), ATOL(1e-8)
      {}
    //! @brief Assignment operator
    template <typename U> Options& operator=
      ( U&options ){
        INTMETH   = options.INTMETH;
        H0        = options.H0;
        HMIN      = options.HMIN;
        HMAX      = options.HMAX;
        RTOL      = options.RTOL;
        ATOL      = options.ATOL;
        return *this;
      }
    //! @brief Enumeration of numerical integration algorithms
    enum METHOD{
      TS=0,	    //!< Explicit Taylor series expansion
      ITS,	    //!< Implicit Taylor series expansion
      RK,		//!< Explicit Runge-Kutta formula
      IRK,		//!< Implicit Runge-Kutta formula
      RAD		//!< Radau collocation
    };
    //! @brief Numerical integration method [Default: MSADAMS]
    METHOD INTMETH;
    //! @brief Initial step-size [Default: 0e0 (auto)]
    double H0;
    //! @brief Minimum step-size [Default: 0e0]
    double HMIN;
    //! @brief Maximum step-size [Default: 0e0 (+inf)]
    double HMAX;
    //! @brief Relative (scalar) integration tolerance
    double RTOL;
    //! @brief Absolute (scalar) integration tolerance
    double ATOL;
  } options;

  //! @brief Structure storing integration statistics
  struct Stats
  {
    //! @brief Constructor
    Stats():
      cputime(0.), numSteps(0), numAE(0)
      {}
    //! @brief Constructor
    void reset()
      { cputime = 0.; numSteps = numAE = 0; }

    //! @brief CPU time
    double cputime;
    //! @brief Integration steps
    unsigned int numSteps;
    //! @brief Calls to algebraic equation solver
    unsigned int numAE;
  };

  // generate <a>N</a> Legendre-Gauss-Radau collocation points between -1 and 1
  static void lgrnodes
    ( unsigned int N, double*x, double eps=1e-7 );
  // generate <a>N</a> Legendre-Gauss-Radau collocation points between <a>t1<\a> and <a>t2</a>
  static void lgrnodes
    ( unsigned int N, double*x, double t1, double t2, double eps=1e-7 );

protected:
  //! @brief Function to initialize GSL statistics
  static void _init_stats
    ( Stats&stats );

  //! @brief Function to finalize GSL statistics
  static void _final_stats
    ( Stats&stats );

  //! @brief Function to display integration statistics
  static void _print_stats
    ( const Stats&stats, std::ostream&os=std::cout );

  //! @brief Private methods to block default compiler methods
  BASE_EXPAND(const BASE_EXPAND&);
  BASE_EXPAND& operator=(const BASE_EXPAND&);
};

inline void
BASE_EXPAND::_init_stats
( Stats&stats )
{
  // Initialize statistics
  stats.reset();
  timeval time;
  gettimeofday(&time, 0) ;
  stats.cputime = - time.tv_sec - time.tv_usec*1e-6;
}

inline void
BASE_EXPAND::_final_stats
( Stats&stats )
{
  // Get final CPU time
  timeval time;
  gettimeofday(&time, 0);
  stats.cputime += time.tv_sec + time.tv_usec*1e-6;
}

inline void
BASE_EXPAND::_print_stats
( const Stats&stats, std::ostream&os )
{
  // Statistics
  os << " No STEPS    " << stats.numSteps
     << std::endl
     << " No AE CALLS " << stats.numAE
     << std::endl
     << " CPU TIME (SEC)     " << std::fixed << std::left
                               << std::setprecision(5) << stats.cputime
     << std::endl << std::endl;
  return;
}

inline void
BASE_EXPAND::lgrnodes
( unsigned int N, double*x, double eps )
{
  if( !N || !x ) return;

  // Use Chebyshev-Gauss-Radau nodes as initial guess for LGR nodes
  for( unsigned int i=0; i<N; i++ )
    x[i] = -std::cos(2*std::acos(-1)*i/(2*(N-1)+1));

  // Iterate correction
  double P[N+1], dP[2], err = eps+1.;
  while( err > eps ){
    err = 0;
    for( unsigned int i=1; i<N; i++ ){
      // 0th and 1st order Legendre polynomials
      P[0] = 1;
      P[1] = x[i];
      // update P from P[2] to P[6]
      for( unsigned int k=1; k<N; k++) // n from 1 to 5
        P[k+1] = ( (2*k+1)*x[i]*P[k] - k*P[k-1] ) / (k+1);
      // calculate derivatives for P[N-1] and P[N]
      dP[1] = N * ( x[i]*P[N] - P[N-1] ) / ( x[i]*x[i] - 1 );       // P'_(N) (x)
      dP[0] = (N-1) * ( x[i]*P[N-1] - P[N-2] ) / ( x[i]*x[i] - 1 ); // P'_(N-1) (x)
      double dx = ( P[N] + P[N-1] ) / ( dP[0] + dP[1] );
      x[i] -= dx; // correction
      if( std::abs(dx) > err ) err = std::abs(dx);
    }
  }
}

inline void
BASE_EXPAND::lgrnodes
( unsigned int N, double*x, double t1, double t2, double eps )
{
  if( !N || !x ) return;
  lgrnodes( N, x, eps );
  for( unsigned int i=0; i<N; i++ ){
    x[i] += 1.; x[i] *= 0.5*(t2-t1);
  }
}

} // end namescape mc

#endif

