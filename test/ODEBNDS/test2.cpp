#include <fstream>
#include "odeslvs_gsl.hpp"
#include "odebnds_gsl.hpp"

#include "interval.hpp"
typedef mc::Interval I;
typedef mc::Ellipsoid E;

#include "tmodel.hpp"
typedef mc::TModel<I> TM;
typedef mc::TVar<I> TV;

int main()
{
      mc::FFGraph IVP;  // DAG describing the problem

      double t0 = 0., tf = 10.;  // Time span
      const unsigned int NS = 50;  // Time stages
      double tk[NS+1]; tk[0] = t0;
      for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + (tf-t0)/(double)NS;

      const unsigned int NP = 2;  // Parameter dimension
      const unsigned int NX = 2;  // State dimension
      const unsigned NF = 2;  // Function dimension

      mc::FFVar P[NP];  // Parameter array
      for( unsigned int i=0; i<NP; i++ ) P[i].set( &IVP );

      mc::FFVar X[NX];  // State array
      for( unsigned int i=0; i<NX; i++ ) X[i].set( &IVP );

      mc::FFVar RHS[NX];  // Right-hand side function
      const double q = 0.1;
      RHS[0] = q*X[0]*(1.-X[0]*X[0]) + X[1]*(1.-q*X[0]*X[1]);
      RHS[1] = q*X[1]*(1.-X[1]*X[1]) - X[0]*(1.+q*X[0]*X[1]) - 0.2*X[1];

      mc::FFVar IC[NX];   // Initial value function
      IC[0] = P[0];
      IC[1] = P[1];

      mc::FFVar FCT[NF];  // 'Objective' function
      FCT[0] = X[0] * X[1];
      FCT[1] = P[0] * pow( X[0], 2 );

      I Ip[NP] = { I(1.75,2.25), I(-0.25,0.25) };
      I *Ixk[NS+1], *Iyk[NS+1];
      for( unsigned k=0; k<=NS; k++ ){
        Ixk[k] = new I[NX];
        Iyk[k] = new I[NF*NX];
      }

      I If[NF], Idf[NF*NP];


      //// SAMPLING
      mc::ODESLVS_GSL<I> LV0;
      LV0.set_dag( &IVP );
      LV0.set_state( NX, X );
      LV0.set_parameter( NP, P );
      LV0.set_differential( NX, RHS );
      LV0.set_initial( NX, IC );
      LV0.set_function( NF, FCT );

      LV0.options.DISPLAY = 1;
      LV0.options.ATOL = LV0.options.RTOL = 1e-8;
      LV0.options.H0 = 1e-6;
      LV0.options.INTMETH = mc::ODESLV_GSL<I>::Options::MSBDF;
      LV0.options.RESRECORD = true;

      // Approximate adjoint bounds
      const unsigned int NSAMP = 20;  // Parameter samples
      for(int g = 0; g < 49; g++){std::cout << "-";}
      std::cout << std::endl;
      std::cout << "|\t\t Approximate Bounds \t\t|" << std::endl; 
      for(int g = 0; g < 49; g++){std::cout << "-";}
      std::cout << std::endl;
      //LV0.bounds( NS, tk, Ip, Ixk, 0, If, NSAMP );
      LV0.bounds_ASA( NS, tk, Ip, Ixk, 0, If, Iyk, Idf, NSAMP );
      std::ofstream apprecSTA("test2_APPROX_STA.dat", std::ios_base::out );
      std::ofstream apprecADJ("test2_APPROX_ADJ.dat", std::ios_base::out );
      LV0.record( apprecSTA, apprecADJ ); 


      //// DIFFERENTIAL INEQUALITIES
      mc::ODEBNDS_GSL<I> LV;
      LV.set_dag( &IVP );
      LV.set_state( NX, X );
      LV.set_parameter( NP, P );
      LV.set_differential( NX, RHS );
      LV.set_initial( NX, IC );
      LV.set_function( NF, FCT );

      LV.options.RESRECORD = true;
      LV.options.DISPLAY = 1;
      LV.options.ATOL = LV.options.RTOL = 1e-8;
      //LV.options.WRAPMIT = mc::ODEBNDS_GSL<I,TM,TV>::Options::DINEQ;
      LV.options.WRAPMIT = mc::ODEBND_GSL<I,TM,TV>::Options::ELLIPS;
      LV.options.HMAX = 1e-2;
      LV.options.H0   = 1e-6;

      // Adjoint bounds from sensitivity class
      for(int g = 0; g < 41; g++){std::cout << "-";}
      std::cout << std::endl;
      std::cout << "|\t State and Adjoint Bounds \t|" << std::endl;     
      for(int g = 0; g < 41; g++){std::cout << "-";}
      std::cout << std::endl;
      LV.bounds_ASA( NS, tk, Ip, Ixk, 0, If, Iyk, Idf );
      std::ofstream direcSTA( "test2_DINEQI_STA.dat", std::ios_base::out );
      std::ofstream direcADJ( "test2_DINEQI_ADJ.dat", std::ios_base::out );
      LV.record( direcSTA, direcADJ );

      // Clean up
      for( unsigned k=0; k<=NS; k++ ){
        delete[] Ixk[k];
        delete[] Iyk[k];
      }

      return 0;
}
