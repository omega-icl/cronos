#include "odeslvs_sundials.hpp"

int main()
{
  mc::FFGraph IVP;  // DAG describing the problem

  double t0 = 0., tf = 10.;   // Time span
  const unsigned int NS = 100;  // Time stages
  double tk[NS+1]; tk[0] = t0;
  for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + (tf-t0)/(double)NS;

  const unsigned NP = 1;  // Parameter dimension
  const unsigned NX = 2;  // State dimension

  mc::FFVar P[NP];  // Parameter array
  for( unsigned i=0; i<NP; i++ ) P[i].set( &IVP );

  mc::FFVar X[NX];  // State array
  for( unsigned i=0; i<NX; i++ ) X[i].set( &IVP );

  mc::FFVar RHS[NX];  // Right-hand side function
  RHS[0] = P[0] * X[0] * ( 1. - X[1] );
  RHS[1] = P[0] * X[1] * ( X[0] - 1. );

  mc::FFVar IC[NX*NS];  // Initial value function
  IC[0] = 1.2;
  IC[1] = 1.1;

  /////////////////////////////////////////////////////////////////////////
  // Compute ODE solutions
  mc::ODESLVS_SUNDIALS LV;

  LV.set_dag( &IVP );
  LV.set_state( NX, X );
  LV.set_parameter( NP, P );
  LV.set_differential( NX, RHS );
  LV.set_initial( NX, IC );

  LV.options.DISPLAY = 1;
  LV.options.RELTOL  = LV.options.ABSTOL = 1e-10;

  double p[NP] = { 3. };  // Parameter values
  double *xk[NS+1];
  for( unsigned k=0; k<=NS; k++ ) xk[k] = new double[NX];
  LV.states( NS, tk, p, xk, 0, 0 );

  for( unsigned k=0; k<=NS; k++ ) delete[] xk[k];

  return 0;
}


