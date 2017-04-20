#include <fstream>
#include "odebnd_sundials.hpp"
#include "interval.hpp"
typedef mc::Interval I;

int main()
{
      mc::FFGraph IVP;  // DAG describing the problem

      const unsigned int NP = 1;  // Parameter dimension
      const unsigned int NX = 2;  // State dimension

      mc::FFVar P[NP];  // Parameter array
      for( unsigned int i=0; i<NP; i++ ) P[i].set( &IVP );

      mc::FFVar X[NX];  // State array
      for( unsigned int i=0; i<NX; i++ ) X[i].set( &IVP );

      mc::FFVar RHS[NX];  // Right-hand side function
      RHS[0] = P[0] * X[0] * ( 1. - X[1] );
      RHS[1] = P[0] * X[1] * ( X[0] - 1. );

      mc::FFVar IC[NX];   // Initial value function
      IC[0] = 1.2;
      IC[1] = 1.1;

      mc::ODEBND_SUNDIALS<I> LV;
      LV.set_dag( &IVP );
      LV.set_time( 0., 25. );
      LV.set_state( NX, X );
      LV.set_parameter( NP, P );
      LV.set_differential( NX, RHS );
      LV.set_initial( NX, IC );

      LV.options.RESRECORD   = true;
      LV.options.ORDMIT      = -2;
      LV.options.DMAX        = 1e2;
      LV.options.NMAX        = 10000;
      LV.options.ODESLV.NMAX = 10000;

      // Compute Interval Bounds
      I Ip[NP];
      Ip[0] = I(2.95,3.05);
      I Ix0[NX], Ixf[NX];
      I* Ixk[2] = {Ix0, Ixf};

      LV.bounds( Ip, Ixk, 0 );
      std::ofstream ofileI( "test0_outer_I.out", std::ios_base::out );
      LV.record( ofileI );

      // Compute Polynomial Bounds
      typedef mc::TModel<I> TM;
      typedef mc::TVar<I>   TV;

      const unsigned int NTM = 5; // Order of Taylor model
      TM TMenv( NP, NTM );        // Taylor model environment
      TV TMp[NP];
      TMp[0] = TV( &TMenv, 0, I(2.95,3.05) );
      TV TMx0[NX], TMxf[NX];
      TV* TMxk[2] = {TMx0, TMxf};

      LV.bounds( TMp, TMxk, 0 );
      std::ofstream ofilePM( "test0_outer_PM.out", std::ios_base::out );
      LV.record( ofilePM );

      // Compute Approximate Bounds (Sampling)
      const unsigned int NSAMP = 50; // Number of sample points
      LV.bounds( Ip, Ixk, 0, NSAMP );
      std::ofstream ofile0( "test0_inner.out", std::ios_base::out );
      LV.record( ofile0 );

      // Compute Hausdorff Distance for Enclosures
      double Hx0[NX], Hxf[NX];
      double* Hxk[2] = {Hx0, Hxf};

      LV.hausdorff( Ip, Hxk, 0, NSAMP );
      LV.hausdorff( TMp, Hxk, 0, NSAMP );

      return 0;
}
