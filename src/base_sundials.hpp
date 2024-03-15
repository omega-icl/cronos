// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__BASE_SUNDIALS_HPP
#define MC__BASE_SUNDIALS_HPP

#include <stdexcept>

#include "sundials/sundials_context.h"
//#include "sundials/sundials_types.h"

namespace mc
{
//! @brief C++ base class for interfacing SUNDIALS solvers
////////////////////////////////////////////////////////////////////////
//! mc::BASE_SUNDIALS is a C++ base class for interfacing SUNDIALS 
//! solvers
////////////////////////////////////////////////////////////////////////
struct BASE_SUNDIALS
{
  //! @brief Default class constructor
  BASE_SUNDIALS
    ()
    { 
      if( SUNContext_Create( SUN_COMM_NULL, &sunctx ) < 0 )
        throw std::runtime_error( "mc::BASE_SUNDIALS: Failed to create SUNContext object" );
    }

  //! @brief Class destructor
  virtual ~BASE_SUNDIALS
    ()
    { 
      SUNContext_Free( &sunctx );  /* Free the SUNDIALS context */
      //if( SUNContext_Free( &sunctx ) < 0 )
      //  throw std::runtime_error( "mc::BASE_SUNDIALS: Failed to free SUNContext object" );
    }

  //! @brief SUNContext object
  SUNContext sunctx;
};

} // end namescape mc

#endif

