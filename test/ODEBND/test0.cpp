#include <fstream>
#include "odeslv_gsl.hpp"
#include "odebnd_gsl.hpp"
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

      mc::ODEBND_GSL<I> LV;
      LV.set_dag( &IVP );
      LV.set_state( NX, X );
      LV.set_parameter( NP, P );
      LV.set_differential( NX, RHS );
      LV.set_initial( NX, IC );
      LV.options.RESRECORD = true;
      LV.options.ORDMIT = 4;

      // Compute Interval Bounds
      I Ip[NP];
      Ip[0] = I(2.95,3.05);
      I Ix0[NX], Ixf[NX];
      I* Ixk[2] = {Ix0, Ixf};
      double tk[2] = {0., 5.};

      LV.bounds( 1, tk, Ip, Ixk );
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

      LV.bounds( 1, tk, TMp, TMxk );
      std::ofstream ofilePM( "test0_outer_PM.out", std::ios_base::out );
      LV.record( ofilePM );

      // Compute Approximate Bounds (Sampling)
       mc::ODESLV_GSL<I> LV0;
      LV0.set_dag( &IVP );
      LV0.set_state( NX, X );
      LV0.set_parameter( NP, P );
      LV0.set_differential( NX, RHS );
      LV0.set_initial( NX, IC );

      const unsigned int NSAMP = 100;  // Number of sample points
      LV0.bounds( 1, tk, Ip, Ixk, 0, 0, NSAMP );
      std::ofstream ofile0( "test0_inner.out", std::ios_base::out );
      LV.record( ofile0 );

      // Compute Hausdorff Distance for Enclosures
      double Hx0[NX], Hxf[NX];
      double* Hxk[2] = {Hx0, Hxf};

      LV.hausdorff( 1, tk, Ip, Hxk, LV0, NSAMP );
      LV.hausdorff( 1, tk, TMp, Hxk, LV0, NSAMP );

      return 0;
}
