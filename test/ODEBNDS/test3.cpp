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

      double t0 = 0., tf = 0.5;  // Time span
      const unsigned int NS = 10;  // Time stages
      double tk[NS+1]; tk[0] = t0;
      for( unsigned k=0; k<NS; k++ ) tk[k+1] = tk[k] + (tf-t0)/(double)NS;

      const unsigned int NP = 1;  // Parameter dimension
      const unsigned int NX = 4;  // State dimension
      const unsigned NF = 2;  // Function dimension

      mc::FFVar P[NP];  // Parameter array
      for( unsigned int i=0; i<NP; i++ ) P[i].set( &IVP );

      mc::FFVar X[NX];  // State array
      for( unsigned int i=0; i<NX; i++ ) X[i].set( &IVP );

      const double g       = 9.81e0;		// m/s^2
      const double m1      = 1.e0;		// kg
      const double m2      = 1.e0;		// kg
      const double l1      = 1.e0;		// m
      const double l2      = 1.e0;		// m
      const double psi1_0  = 3.*mc::PI/4.;	// rad
      const double psi2_0  = -11.*mc::PI/20.;	// rad
      const double psi3_0  = .43e0;		// rad/s
      const double psi4_0  = .67e0;		// rad/s

      mc::FFVar IC[NX];  // Initial value function
      IC[0] = psi1_0 * P[0];
      IC[1] = psi2_0;
      IC[2] = psi3_0;
      IC[3] = psi4_0;

      mc::FFVar RHS[NX];  // RHS of differential equations
      mc::FFVar a11 = m1*l1 + m2*(l1+l2*cos(X[1]));
      mc::FFVar a21 = m2*(l1*cos(X[1])+l2);
      mc::FFVar a12 = m2*l2*cos(X[1]);
      mc::FFVar a22 = m2*l2;
      mc::FFVar b1  = -g*(m1+m2)*sin(X[0]) + m2*l2*sin(X[1])*((X[2]+X[3])*(X[2]+X[3]));
      mc::FFVar b2  = -g*m2*sin(X[0]+X[1]) - m2*l1*sin(X[1])*(X[2]*X[2]);
      RHS[0] = X[2];
      RHS[1] = X[3];
      RHS[2] = ( a22*b1 - a12*b2 ) / ( a11*a22 - a12*a21 );
      RHS[3] = ( a11*b2 - a21*b1 ) / ( a11*a22 - a12*a21 );

      mc::FFVar FCT[NF];  // 'Objective' function
      FCT[0] = X[0] * X[1];
      FCT[1] = P[0] * X[3] * X[4];

      I Ip[NP] = { I(0.995,1.005) };
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
      LV0.options.ATOL = LV0.options.RTOL = 1e-10;
      LV0.options.INTMETH = mc::ODESLV_GSL<I>::Options::MSBDF;
      LV0.options.RESRECORD = true;

      // Approximate adjoint bounds
      const unsigned int NSAMP = 20;  // Parameter dimension
      for(int g = 0; g < 45; g++){std::cout << "-";}
      std::cout << std::endl;
      std::cout << "|\t\t Sampled Bounds \t\t|" << std::endl; 
      for(int g = 0; g < 45; g++){std::cout << "-";}
      std::cout << std::endl;
      //LV0.bounds( NS, tk, Ip, Ixk, 0, If, NSAMP );
      LV0.bounds_ASA( NS, tk, Ip, Ixk, 0, If, Iyk, Idf, NSAMP );
      std::ofstream apprecSTA("test3_APPROX_STA.dat", std::ios_base::out );
      std::ofstream apprecADJ("test3_APPROX_ADJ.dat", std::ios_base::out );
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

      // Adjoint bounds from sensitivity class
      for(int g = 0; g < 41; g++){std::cout << "-";}
      std::cout << std::endl;
      std::cout << "|\t State and Adjoint Bounds \t|" << std::endl;     
      for(int g = 0; g < 41; g++){std::cout << "-";}
      std::cout << std::endl;
      LV.bounds_ASA( NS, tk, Ip, Ixk, 0, If, Iyk, Idf );
      std::ofstream direcSTA( "test3_DINEQI_STA.dat", std::ios_base::out );
      std::ofstream direcADJ( "test3_DINEQI_ADJ.dat", std::ios_base::out );
      LV.record( direcSTA, direcADJ );

      // Clean up
      for( unsigned k=0; k<=NS; k++ ){
        delete[] Ixk[k];
        delete[] Iyk[k];
      }

      return 0;
}
