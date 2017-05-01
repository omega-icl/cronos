#include <fstream>
#include "odebnd_sundials.hpp"
#include "interval.hpp"
typedef mc::Interval I;
typedef mc::CModel<I> TM;
typedef mc::CVar<I>   TV;

int main()
{
  mc::FFGraph IVP;  // DAG describing the problem

  const unsigned int NP = 2;  // Parameter dimension
  const unsigned int NX = 2;  // State dimension

  mc::FFVar P[NP];  // Parameter array
  for( unsigned int i=0; i<NP; i++ ) P[i].set( &IVP );

  mc::FFVar X[NX];  // State array
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &IVP );

  mc::FFVar RHS[NX];  // Right-hand side function
  RHS[0] = P[0] * X[0] * ( 1. - X[1] );
  RHS[1] = P[1] * X[1] * ( X[0] - 1. );

  mc::FFVar IC[NX];   // Initial value function
  IC[0] = 1.2;
  IC[1] = 1.1;

  mc::ODEBND_SUNDIALS<I,TM,TV> LV;
  LV.set_dag( &IVP );
  LV.set_time( 0., 60. );
  LV.set_state( NX, X );
  LV.set_parameter( NP, P );
  LV.set_differential( NX, RHS );
  LV.set_initial( NX, IC );

  LV.options.RESRECORD   = 1000;
  LV.options.ORDMIT  = -2;
  LV.options.DMAX    = 1e2;
  LV.options.NMAX    = 10000;
  LV.options.ODESLVS     = LV.options;

  // Compute Interval Bounds
  I Ip[NP];
  Ip[0] = I(2.99,3.01);
  Ip[1] = I(0.99,1.01);

  LV.bounds( Ip );
  std::ofstream ofileI( "test0_DINEQI_STA.dat", std::ios_base::out );
  LV.record( ofileI );

  // Compute Polynomial Bounds
  const unsigned int NTM = 5; // Order of Taylor model
  TM TMenv( NP, NTM );    // Taylor model environment
  TV TMp[NP];
  TMp[0] = TV( &TMenv, 0, I(2.99,3.01) );
  TMp[1] = TV( &TMenv, 1, I(0.99,1.01) );

  LV.bounds( TMp );
  std::ofstream ofilePM( "test0_DINEQPM_STA.dat", std::ios_base::out );
  LV.record( ofilePM );

  // Compute Approximate Bounds (Sampling)
  const unsigned int NSAMP = 50; // Number of sample points
  LV.bounds( NSAMP, Ip );
  std::ofstream ofile0( "test0_APPROX_STA.dat", std::ios_base::out );
  LV.record( ofile0 );

  // Compute Hausdorff Distance for Enclosures
  double Hx0[NX], Hxf[NX];
  double* Hxk[2] = {Hx0, Hxf};

  LV.set_time( 0., 8. );
  LV.hausdorff( NSAMP, Ip, Hxk );
  LV.hausdorff( NSAMP, TMp, Hxk );

  return 0;
}
