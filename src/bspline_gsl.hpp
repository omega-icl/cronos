// Copyright (C) 2016 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__BSPLINE_GSL_HPP
#define MC__BSPLINE_GSL_HPP

#include <iostream>
#include <vector>

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_statistics.h>

namespace mc
{
//! @brief C++ class for B-spline approximation of sampled data using GSL
////////////////////////////////////////////////////////////////////////
//! mc::BSPLINE_GSL is a C++ class for B-spline approximation of sampled
//! data using GSL.
////////////////////////////////////////////////////////////////////////
class BSPLINE_GSL
{
public:
  //! @brief sampled data 
  std::vector< std::vector<double> > data;

private:
  //! @brief GSL workspace for B-spline basis
  gsl_bspline_workspace* _bw;
  //! @brief GSL workspace for B-spline fits
  std::vector<gsl_multifit_linear_workspace*> _mw;
  //! @brief Number of breakpoints
  unsigned _nk;
  //! @brief Order of B-splines
  unsigned _nord;
  //! @brief B-spline basis
  gsl_vector* _B;
  //! @brief Fit matrix
  gsl_matrix* _X;
  //! @brief Data vector
  gsl_vector* _y;
  //! @brief Fit result
  std::vector<gsl_vector*> _c;
  //! @brief Covariance matrix
  std::vector<gsl_matrix*> _V;

public:
  //! @brief Constructor
  BSPLINE_GSL
    ();
  //! @brief Destructor
  ~BSPLINE_GSL
    ();
  //! @brief Initialize B-spline with <a>ndim</a> dimensions
  bool set
    ( const unsigned NDIM, const unsigned PREALLOC, const unsigned NK,
      const unsigned NORD=4 );
  //! @brief Add new data (must be of dimension <a>ndim</a>)
  bool add
    ( const double&t, const double*x );
  //! @brief Fit B-spline to data
  bool fit
    ();
  //! @brief Evaluate data (all <a>ndim</a> dimensions) at time <a>t</a> in stage <a>istg</a> of mesh
  bool eval
    ( const double t, double*x );
};

inline
BSPLINE_GSL::BSPLINE_GSL
()
: _bw(0), _nk(0), _nord(0), _B(0), _X(0), _y(0)
{}

inline
BSPLINE_GSL::~BSPLINE_GSL
()
{
  gsl_vector_free( _B );
  gsl_vector_free( _y );
  gsl_matrix_free( _X );
  for( auto it=_c.begin(); it!=_c.end(); ++it )
    gsl_vector_free( *it );
  for( auto it=_V.begin(); it!=_V.end(); ++it )
    gsl_matrix_free( *it );
  gsl_bspline_free( _bw );
  for( auto it=_mw.begin(); it!=_mw.end(); ++it )
    gsl_multifit_linear_free( *it );
}

inline bool
BSPLINE_GSL::set
( const unsigned NDIM, const unsigned PREALLOC, const unsigned NK,
  const unsigned NORD )
{
  try{
    data.assign( NDIM+1, std::vector<double>() );
    for( auto it=data.begin(); it!=data.end(); ++it )
      it->reserve( PREALLOC );

    _nk = NK; _nord = NORD;
    gsl_bspline_free( _bw );
    _bw = gsl_bspline_alloc( NORD, NK );
  }
  catch(...){
    return false;
  }
  return true;
}

inline bool
BSPLINE_GSL::add
( const double&t, const double*x )
{
  try{
    data[0].push_back( t );
    for( unsigned i=0; i<data.size()-1; i++ )
      data[i+1].push_back( x[i] );
  }
  catch(...){
    return false;
  }
  return true;
}

inline bool
BSPLINE_GSL::fit
()
{
  try{
    for( auto it=_mw.begin(); it!=_mw.end(); ++it )
      gsl_multifit_linear_free( *it );
    _mw.clear();

    // initialize breakpoints and linear fit
    const unsigned NDIM = data.size()-1;
    const unsigned NP   = data[0].size();
    const unsigned NC   = _nk+_nord-2;
    for( unsigned i=0; i<NDIM; i++ )
      _mw.push_back( gsl_multifit_linear_alloc( NP, NC ) );

    // construct fit matrix X
    gsl_bspline_knots_uniform( data[0][0], data[0][NP-1], _bw );
    gsl_vector_free(_B); _B = gsl_vector_alloc( NC );
    gsl_matrix_free(_X); _X = gsl_matrix_alloc( NP, NC );
    for( unsigned k=0; k<NP; ++k ){
      gsl_bspline_eval( data[0][k], _B, _bw);
      for( unsigned j=0; j<NC; ++j )
        gsl_matrix_set( _X, k, j, gsl_vector_get( _B, j ) );
    }

    // do the fit
    double ssr;
    gsl_vector_free( _y ); _y = gsl_vector_alloc( NP );
    for( auto it=_c.begin(); it!=_c.end(); ++it )
      gsl_vector_free( *it );
    _c.clear();
    for( auto it=_V.begin(); it!=_V.end(); ++it )
      gsl_matrix_free( *it );
    _V.clear();
    for( unsigned i=0; i<NDIM; i++ ){
      _c.push_back( gsl_vector_alloc( NC ) );
      _V.push_back( gsl_matrix_alloc( NC, NC ) );
      for( unsigned k=0; k<NP; ++k )
        gsl_vector_set( _y, k, data[i+1][k] );
      gsl_multifit_linear( _X, _y, _c[i], _V[i], &ssr, _mw[i] );
      std::cerr << "SSR" << i << " = " << ssr << std::endl;
    }
  }
  catch(...){
    return false;
  }
  return true;
}

inline bool
BSPLINE_GSL::eval
( const double t, double*x )
{
  try{
    const unsigned NDIM = data.size()-1;
    double xstd;
    gsl_bspline_eval( t, _B, _bw );
    for( unsigned i=0; i<NDIM; i++ )
      gsl_multifit_linear_est( _B, _c[i], _V[i], x+i, &xstd );
  }
  catch(...){
    return false;
  }
  return true;
}

} // mc namespace
#endif

