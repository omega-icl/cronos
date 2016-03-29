// Copyright (C) 2015 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__MESH_GSL_HPP
#define MC__MESH_GSL_HPP

#include <iostream>
#include <vector>

#include "gsl/gsl_errno.h"
#include "gsl/gsl_interp.h"

namespace mc
{
//! @brief C++ class storing and processing data mesh using GSL.
////////////////////////////////////////////////////////////////////////
//! mc::MESH_GSL is a C++ class for storing and processing data mesh
//! using the interpolation method gsl_interp in GSL.
////////////////////////////////////////////////////////////////////////
class MESH_GSL
{
public:
  //! @brief Enumeration of interpolation algorithms
  enum INTERPOLATION_METHOD{
    LINEAR=0,	//!< Linear interpolation. [Default]
    CSPLINE,	//!< Cubic spline with natural boundary conditions.
    AKIMA	//!< Non-rounded Akima spline with natural boundary conditions.
  };
  //! @brief state bound parameterizations, one for every ODE and time point 
  std::vector< std::vector<double> > data;
  //! @brief starting index of every stage in <a>data</a>
  std::vector<unsigned> index;
  //! @brief type of state bound parameterization in <a>data</a>
  short type;

private:
  //! @brief GSL lookup accelerators for interpolation
  std::vector<gsl_interp_accel*> _accel;
  //! @brief GSL drivers for interpolation
  std::vector<gsl_interp*> _driver;
  //! @brief Resize driver
  void _resize
    ( const int driver_type, const unsigned ndat );

public:
  //! @brief Destructor
  ~MESH_GSL();
  //! @brief Initialize mesh with <a>ndim</a> dimensions over <a>nstg</a> stages
  bool set
    ( const unsigned nstg, const unsigned ndim, const unsigned prealloc,
      const short info=0 );
  //! @brief Add new data (must be of dimension <a>ndim</a>) to stage <a>istg</a> of mesh
  bool add
    ( const unsigned istg, const double&t, const double*x, const bool index=false );
  //! @brief Interpolate data (all <a>ndim</a> dimensions) in stage <a>istg</a> of mesh
  bool interp
    ( const unsigned istg, const INTERPOLATION_METHOD type=CSPLINE );
  //! @brief Evaluate data (all <a>ndim</a> dimensions) at time <a>t</a> in stage <a>istg</a> of mesh
  bool eval
    ( const unsigned istg, const double t, double*x );
};

inline
MESH_GSL::~MESH_GSL
()
{
  for( unsigned i=0; i<_driver.size(); i++ )
    gsl_interp_free( _driver[i] );
  for( unsigned i=0; i<_accel.size(); i++ )
    gsl_interp_accel_free( _accel[i] );
}

inline bool
MESH_GSL::set
( const unsigned nstg, const unsigned ndim, const unsigned prealloc,
  const short info )
{
  try{
    // (Re)initialize mesh data and indexing
    data.assign( ndim+1, std::vector<double>() );
    for( auto it=data.begin(); it!=data.end(); ++it ) it->reserve(prealloc);
    index.assign( nstg, 0 );
  }
  catch(...){
    return false;
  }
  type = info;
  return true;
}

inline bool
MESH_GSL::add
( const unsigned istg, const double&t, const double*x, const bool addindex )
{
  if( index.size()<istg+1 ) return false;

  // Index current mesh position
  if( addindex ) index[istg] = data[0].size();

  // Store time and state into mesh
  data[0].push_back( t );
  for( unsigned i=0; i<data.size()-1; i++ ) data[i+1].push_back( x[i] );
  return true;
}

inline bool
MESH_GSL::interp
( const unsigned istg, const INTERPOLATION_METHOD type )
{
  const unsigned offset = index[istg-1];
  const unsigned ndat = ( istg==index.size()? data[0].size(): index[istg] ) - offset; 

  // Initialize interpolation driver and accelerator
  _resize( type, ndat );

  // Interpolate every state in current stage
  for( unsigned i=0; i<_driver.size(); i++ ){
#ifdef MC__ODEBND_GSL_INTERP_DEBUG
    std::cout << "Stage #" << istg << "  State #" << i << std::endl; 
    std::cout << data[0].data() << "  " << data[i+1].data() << std::endl;
    for( unsigned k=0; k<ndat; k++ )
      std::cout << (data[0].data()+offset)[k] << "  "
                << (data[i].data()+offset)[k] << std::endl;
    { int dum; std::cin >> dum; }
#endif
    int flag = gsl_interp_init( _driver[i], data[0].data()+offset, data[i].data()+offset, ndat );
#ifdef MC__ODEBND_GSL_INTERP_DEBUG
    std::cout << "flag: " << flag << std::endl;
#endif
    if( flag != GSL_SUCCESS ) return false;
  }
  // !!!gsl_interp_init returns an int flag, but it is unclear what it is!!!
  // Function: int gsl_interp_init (gsl_interp * interp, const double xa[], const double ya[], size_t size)

  // Reset lookup accelerator
  if( _accel.size() < data.size()-1 ) return false;
  for( unsigned i=0; i<_accel.size(); i++ )
    if( gsl_interp_accel_reset( _accel[i] ) != GSL_SUCCESS ) return false;
  return true;
}

inline void
MESH_GSL::_resize
( const int driver_type, const unsigned ndat  )
{
  // Set GSL interpolation lookup accelerator - resize only if necessary
  const unsigned naccel = _accel.size();
  for( unsigned i=data.size()-1; i<naccel; i++ ) gsl_interp_accel_free( _accel[i] );
  _accel.resize( data.size()-1 );
  for( unsigned i=naccel; i<data.size()-1; i++ ) _accel[i] = gsl_interp_accel_alloc();

  // Reset GSL interpolation drivers
  for( unsigned i=0; i<_driver.size(); i++ ) gsl_interp_free( _driver[i] );
  _driver.clear();
  switch( driver_type ){
  case AKIMA:
    if( ndat >= 5 ){
      for( unsigned i=0; i<data.size()-1; i++ ) _driver.push_back( gsl_interp_alloc( gsl_interp_akima, ndat ) );
      break;
    }
  case CSPLINE:
    if( ndat >= 3 ){
      for( unsigned i=0; i<data.size()-1; i++ ) _driver.push_back( gsl_interp_alloc( gsl_interp_cspline, ndat ) );
      break;
    }
  case LINEAR: default:
    for( unsigned i=0; i<data.size()-1; i++ ) _driver.push_back( gsl_interp_alloc( gsl_interp_linear, ndat ) );
    break;
  }
}

inline bool
MESH_GSL::eval
( const unsigned istg, const double t, double*x )
{
  if( _driver.size() != data.size()-1 ) return false;
  const unsigned offset = index[istg-1];
  for( unsigned i=0; i<_driver.size(); i++ )
    if( gsl_interp_eval_e( _driver[i], data[0].data()+offset, data[i+1].data()+offset, t, _accel[i], x+i ) == GSL_EDOM )
      return false;
  return true;
}

} // mc namespace
#endif
