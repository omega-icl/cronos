#include <fstream>
#include "odebnd_expand.hpp"
#include "interval.hpp"
#include "scmodel.hpp"

typedef mc::Interval I;
//typedef mc::SCModel<I> PMI;
//typedef mc::SCVar<I>   PVI;
typedef mc::CModel<I> PMI;
typedef mc::CVar<I>   PVI;

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
      RHS[0] = ( 3. + P[0] ) * X[0] * ( 1. - X[1] );
      RHS[1] = ( 3. + P[0] ) * X[1] * ( X[0] - 1. );

      mc::FFVar IC[NX];   // Initial value function
      IC[0] = 1.2;// + 0.1*P[0];
      IC[1] = 1.1;// + P[0];

      mc::ODEBND_EXPAND<I,PMI,PVI> LV;
      LV.options.H0 = 0.1;
      LV.options.LBLK = 20;
      LV.options.DBLK = 20;
      LV.options.RESRECORD = true;
      LV.options.DISPLAY = 2;

      LV.options.AEBND.MAXIT = 100;
      LV.options.AEBND.DISPLAY = 1;
      LV.options.AEBND.BLKDEC  = false;
      LV.options.AEBND.RTOL    =
      LV.options.AEBND.ATOL    = 1e-10;
      LV.options.AEBND.BOUNDER = mc::AEBND<I,PMI,PVI>::Options::GS;//KRAW;//AUTO;
      LV.options.AEBND.PRECOND = mc::AEBND<I,PMI,PVI>::Options::INVMB;//INVMB;//INVBD;//NONE;

      LV.set_dag( &IVP );
      LV.set_state( NX, X );
      LV.set_parameter( NP, P );
      LV.set_differential( NX, RHS );
      LV.set_initial( NX, IC );
      LV.setup( 1 );

      // Compute Interval Bounds
      I Ip[NP];
      Ip[0] = 0.05 * I(-1.,1);
      I Ix0[NX], Ixf[NX];
      I* Ixk[2] = {Ix0, Ixf};
      double tk[2] = {0., 25. };

      //LV.bounds( 1, tk, Ip, Ixk, 0 );
      std::ofstream ofileI( "test0_EXPANDI_STA.dat", std::ios_base::out );
      LV.record( ofileI );

      // Compute Polynomial Bounds
      const unsigned int NPM = 5; // Order of polynomial model
      PMI PMenv( NP, NPM );       // Polynomial model environment
      //PMI PMenv( NPM );       // Polynomial model environment
      PVI PMp[NP];
      PMp[0] = PVI( &PMenv, 0, Ip[0] );
      //PVI PMx0[NX] = { I(0.5,1.5), I(0.5,1.5) }, PMxf[NX] = { I(0.5,1.5), I(0.5,1.5) };
      PVI PMx0[NX] = { I(0.,5.), I(0.,5.) }, PMxf[NX] = { I(0.,5.), I(0.,5.) };
      PVI* PMxk[2] = {PMx0, PMxf};

      LV.bounds( 1, tk, PMp, PMxk, 0 );
      std::ofstream ofilePM( "test0_EXPANDPM_STA.dat", std::ios_base::out );
      LV.record( ofilePM );

      // Compute Approximate Bounds (Sampling)
      const unsigned int NSAMP = 50; // Number of sample points
      LV.bounds( 1, tk, Ip, Ixk, 0, NSAMP );
      std::ofstream ofile0( "test0_APPROX_STA.dat", std::ios_base::out );
      LV.record( ofile0 );

      // Compute Hausdorff Distance for Enclosures
      double Hx0[NX], Hxf[NX];
      double* Hxk[2] = {Hx0, Hxf};

      //LV.hausdorff( 1, tk, Ip, Hxk, 0, NSAMP );
      //LV.hausdorff( 1, tk, TMp, Hxk, 0, NSAMP );

      return 0;
}
