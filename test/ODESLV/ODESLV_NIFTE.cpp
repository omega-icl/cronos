// ODESLV_NIFTE.cpp -- repurposed from test3.cpp (NIFTE dynamics, NTP model).
// =============================================================================================
//   5 states (P_ad, P_d, P_p, U_f, U_p), 10 design parameters (R_f .. C_p), 5 constants (P0, U0, tau, Lambda, K),
//   t in [0,5], one stage.  Functions at tf: U_l and U_th (the quadrature integrands evaluated at tf, as test3).
//   The parameters span 1e-9 .. 5e8, so gradients are compared as SCALED sensitivities p_i dF_k/dp_i, relative to
//   the largest one of each function.  The model is DECLARED through FFModel (add_state / add_input / set_constant /
//   add_equation / add_output); test3's quadrature states are not used by its functions and are dropped.
//
// PART A -- ODESLVS on the declared model, SPARSE linear solver, values given by name.
//     A1  solve_state, solve_sensitivity, solve_adjoint all NORMAL
//     A2  forward == adjoint over all 10 parameters (threshold 1e-4, set after measuring; LV: 1.5e-7)
// PART B -- FFODESLV on the SAME declared model (two-map form):
//   map 1 = the resistances { R_f, R_l, R_th };  map 2 = the other 7 parameters and the 5 constants.
//     B1  op value == solve_state
//     B2  numeric (FFGradODESLV) dF/d(map 1) == forward sensitivity;  map-2 columns ZERO by contract
//     B3  SYMDIFF = { K (constant), L_d (map-2 parameter) }: dF/dL_d == forward sensitivity,
//         dF/dK == central FD (named solves at 1e-13);  every other column ZERO by contract
//     B4  SYMDIFF = map 1, the SECOND SYMDIFF on the same op: == forward sensitivity (scaled)
//
// Build: g++ -std=c++17 ... ODESLV_NIFTE.cpp

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include "ffode.hpp"

static int npass = 0, nfail = 0;
static void check( bool c, std::string const& w ){ std::cout << "  " << ( c? "PASS  ": "FAIL  " ) << w << "\n"; c? ++npass: ++nfail; }

double const t0 = 0., tf = 5.;
std::vector<double> const cval{ 1.013e5, 8.00e-4, 5., 330., 1.37 };                                  // P0 U0 tau Lambda K
std::vector<double> const pval{ 2.13e6, 4.08e6, 5.02e8, 1.80e5, 4.74e6, 1.27e7, 3.77e5, 1.76e-9, 7.35e-8, 7.43e-8 };
std::vector<double> const xic { 0.01, 0., 0.01, 0., 0. };
char const* const PN[10] = { "R_f", "R_l", "R_th", "L_d", "L_f", "L_l", "L_p", "C_ad", "C_d", "C_p" };
char const* const CN[5]  = { "P0", "U0", "tau", "Lambda", "K" };
char const* const XN[5]  = { "P_ad", "P_d", "P_p", "U_f", "U_p" };
enum { R_f, R_l, R_th, L_d, L_f, L_l, L_p, C_ad, C_d, C_p };
enum { P0, U0, TAU, LAMBDA, K };

//! @brief NIFTE right-hand sides, quadrature integrands and the two functions U_l, U_th (test3's expressions).
static void nifte( std::vector<mc::FFVar> const& X, std::vector<mc::FFVar> const& P, std::vector<mc::FFVar> const& C,
                   std::vector<mc::FFVar>& RHS, std::vector<mc::FFVar>& QUAD, std::vector<mc::FFVar>& F )
{
  auto const &Pad = X[0], &Pd = X[1], &Pp = X[2], &Uf = X[3], &Up = X[4];
  auto const &Rf = P[R_f], &Rl = P[R_l], &Rth = P[R_th], &Ld = P[L_d], &Lf = P[L_f], &Ll = P[L_l], &Lp = P[L_p],
             &Cad = P[C_ad], &Cd = P[C_d], &Cp = P[C_p];
  auto const &p0 = C[P0], &u0 = C[U0], &tau = C[TAU], &Lam = C[LAMBDA], &Kg = C[K];
  mc::FFVar const den = Ll*(Ld+Lf) + Lp*(Ld+Lf+Ll);
  RHS = { ( tau*(Kg*tanh(Lam*Pd)-Pad) )/(Rth*Cad) + ( tau*u0*(Uf + Up) )/(p0*Cad),
          ( tau*u0*Uf )/(p0*Cd),
          ( tau*u0*Up )/(p0*Cp),
          ( (tau*p0/u0)*( Ll*Pp-Lp*Pad-(Lp+Ll)*Pd ) - tau*(Ll*Rf+Lp*(Rf+Rl))*Uf - tau*Lp*Rl*Up ) / den,
          ( (tau*p0/u0)*( Ll*Pd - (Ld+Lf)*Pad - (Ld+Lf+Ll)*Pp ) + tau*(Ll*Rf-(Ld+Lf)*Rl)*Uf - tau*(Ld+Lf)*Rl*Up ) / den };
  mc::FFVar const dPaddt = ( tau*(Kg*((mc::exp(2*Lam*Pd)-1.)/(mc::exp(2*Lam*Pd)+1.))-Pad) )/(Rth*Cad) + ( tau*u0*(Uf + Up) )/(p0*Cad);
  mc::FFVar const Uad = Cad*dPaddt*p0/tau;
  QUAD = { Uf*u0, -(Uf + Up)*u0, Uad - Uf*u0 - Up*u0 };
  F    = { QUAD[1], QUAD[2] };                                                        // U_l, U_th (as test3)
}

//! @brief worst over functions k of max_i |a_ik - b_ik| |p_i| / max_i |b_ik p_i|   (G[i][k], i = parameter)
static double scaled_err( std::vector<std::vector<double>> const& a, std::vector<std::vector<double>> const& b,
                          std::vector<double> const& p )
{
  double w = 0.;
  for( size_t k = 0; k < b[0].size(); ++k ){
    double s = 0.;  for( size_t i = 0; i < p.size(); ++i ) s = std::max( s, std::fabs( b[i][k] * p[i] ) );
    for( size_t i = 0; i < p.size(); ++i ) w = std::max( w, std::fabs( a[i][k] - b[i][k] ) * std::fabs( p[i] ) / s );
  }
  return w;
}
static void tol_opts( mc::ODESLVS_CVODES& I, double tol, bool sparse )
{
  I.options.INTMETH = mc::BASE_CVODES::Options::MSBDF;  I.options.NLINSOL = mc::BASE_CVODES::Options::NEWTON;
  I.options.LINSOL  = sparse? mc::BASE_CVODES::Options::SPARSE: mc::BASE_CVODES::Options::DENSE;
  I.options.DISPLAY = 0;
  // NMAX: at RTOL 1e-10 the stiff stretch around t = 3.5 needs more than test3's 10000 steps (mxstep otherwise).
  I.options.NMAX = 1000000;
  I.options.ATOL = I.options.ATOLB = I.options.ATOLS = tol;  I.options.RTOL = I.options.RTOLB = I.options.RTOLS = tol;
}

int main()
{
  std::cout << std::scientific << std::setprecision(6)
            << "================================================================\n"
            << "  ODESLV_NIFTE: NIFTE dynamics -- ODESLVS and FFODESLV on one model declared through FFModel\n"
            << "================================================================\n";
  size_t const np = 10, nf = 2;

  mc::FFGraph D;  mc::ODESLVS_CVODES I( &D );
  mc::FFVar t = D.add_var( "t" );
  std::vector<mc::FFVar> Xd( 5 ), Pd( np ), Cd( 5 );
  for( size_t i = 0; i < 5; ++i )  Xd[i] = D.add_var( std::string( XN[i] ) + "(t)" );
  for( size_t i = 0; i < np; ++i ) Pd[i] = D.add_var( PN[i] );
  for( size_t i = 0; i < 5; ++i )  Cd[i] = D.add_var( CN[i] );
  std::vector<mc::FFVar> RHSd, QUADd, Fd;  nifte( Xd, Pd, Cd, RHSd, QUADd, Fd );
  mc::FFPartial OpP;  mc::FFEval OpE;
  I.add_domain( t, mc::FFDom( t0, tf, 1, mc::FFDom::LGR, 4 ) );  I.set_evolution_domain( t );
  for( size_t i = 0; i < 5; ++i ){ I.add_state( Xd[i], {t} ); I.update_ref( Xd[i], xic[i] ); }
  for( size_t i = 0; i < np; ++i ) I.add_input( Pd[i] );
  I.FFModel::set_constant( Cd, cval );
  int const T_INT = mc::FFDom::ALL - mc::FFDom::LB;
  auto const io = mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INTERIOR, 0 );
  auto const ii = mc::FFModel::EqnOptions( mc::FFModel::EqnRole::INITIAL,  0 );
  for( size_t i = 0; i < 5; ++i ){
    I.add_equation( OpP( Xd[i], t ) - RHSd[i], {t}, {T_INT},         io );
    I.add_equation( Xd[i] - xic[i],             {t}, {mc::FFDom::LB}, ii );
  }
  for( auto const& f : Fd ) I.add_output( OpE( f, t, tf ) );              // U_l, U_th at tf (as test3)
  tol_opts( I, 1e-10, true );
  std::vector<mc::FFModel::InputVal> vIn;                                    // every input and constant, by name
  for( size_t i = 0; i < np; ++i ) vIn.emplace_back( Pd[i], std::vector<double>{ pval[i] } );
  for( size_t i = 0; i < 5; ++i )  vIn.emplace_back( Cd[i], std::vector<double>{ cval[i] } );

  std::cout << "\n--- PART A: ODESLVS on the declared model\n";
  check( I.setup(), "A0 declared model sets up" );
  auto const s0 = I.solve( vIn );         auto const FB  = I.val_function();
  auto const s1 = I.solve_fsens( vIn );   auto const GfA = I.val_function_gradient();
  auto const s2 = I.solve_asens( vIn );       auto const GaA = I.val_function_gradient();
  bool const okA = s0 == mc::ODESLVS_CVODES::STATUS::NORMAL && s1 == mc::ODESLVS_CVODES::STATUS::NORMAL && s2 == mc::ODESLVS_CVODES::STATUS::NORMAL;
  std::cout << "  F = ( U_l, U_th )(tf) = ( " << FB[0] << ", " << FB[1] << " )\n";
  check( okA && GfA.size() == np && GaA.size() == np, "A1 solve_state, solve_sensitivity, solve_adjoint NORMAL; 10 gradient rows" );
  if( !okA ){ std::cout << "\n  " << npass << " passed, " << nfail << " failed -- FAILURES\n"; return 1; }
  { std::vector<double> pord( np );                                        // parameter values in var_parameter() order
    for( size_t i = 0; i < np; ++i ) pord[ I.parameter_index( Pd[i] )[0] ] = pval[i];
    std::ostringstream o; double const w = scaled_err( GaA, GfA, pord );
    o << "A2 forward == adjoint, 10 parameters (worst scaled " << w << ")"; check( w < 1e-4, o.str() ); }

  std::cout << "\n--- PART B: FFODESLV on the same declared model\n";
  I.solve_fsens( vIn );  auto const Gref = I.val_function_gradient();  // rows by parameter index
  auto row = [&]( size_t i ){ return Gref[ I.parameter_index( Pd[i] )[0] ]; };

  mc::FFGraph NLP;
  std::vector<mc::FFVar> pv( np ), cv( 5 );
  for( size_t i = 0; i < np; ++i ) pv[i] = NLP.add_var( std::string( "v" ) + PN[i] );
  for( size_t i = 0; i < 5; ++i )  cv[i] = NLP.add_var( std::string( "v" ) + CN[i] );
  std::vector<mc::FFModel::InputArg> map1, map2;
  for( size_t i : { R_f, R_l, R_th } ) map1.emplace_back( Pd[i], std::vector<mc::FFVar>{ pv[i] } );
  for( size_t i = 0; i < np; ++i ) if( i != R_f && i != R_l && i != R_th ) map2.emplace_back( Pd[i], std::vector<mc::FFVar>{ pv[i] } );
  for( size_t i = 0; i < 5; ++i ) map2.emplace_back( Cd[i], std::vector<mc::FFVar>{ cv[i] } );
  mc::FFODESLV op;  std::vector<mc::FFVar> Fop;
  try{ Fop = op( map1, map2, &I, mc::FFODESLV::COPY, "NIFTE" ); }
  catch( std::exception& e ){ check( false, std::string( "B1 two-map embedding threw: " ) + e.what() ); return 1; }
  std::vector<mc::FFVar> vX( pv );  vX.insert( vX.end(), cv.begin(), cv.end() );   // 15 columns: 10 params, 5 constants
  std::vector<double> vXv( pval );  vXv.insert( vXv.end(), cval.begin(), cval.end() );
  size_t const nX = vX.size();
  std::vector<double> Fv( nf );  NLP.eval( Fop, Fv, vX, vXv );
  { double w = 0.; for( size_t k = 0; k < nf; ++k ) w = std::max( w, std::fabs( Fv[k] - FB[k] ) / std::fabs( FB[k] ) );
    std::ostringstream o; o << "B1 op value == solve_state (worst rel " << w << ")"; check( w < 1e-9, o.str() ); }

  auto jac = [&]( std::vector<mc::FFVar> const& sym ){                         // d[k*nX + col]
    mc::FFODESLV::options.SYMDIFF = sym;
    auto const dF = NLP.FAD( Fop, vX );  std::vector<double> d( dF.size() );  NLP.eval( dF, d, vX, vXv );
    mc::FFODESLV::options.SYMDIFF.clear();  return d; };
  // scaled error of the op's column for parameter i against the reference row i
  auto col_err = [&]( std::vector<double> const& d, std::vector<size_t> const& cols ){
    double w = 0.;
    for( size_t k = 0; k < nf; ++k ){
      double s = 0.;  for( size_t i = 0; i < np; ++i ) s = std::max( s, std::fabs( row( i )[k] * pval[i] ) );
      for( size_t i : cols ) w = std::max( w, std::fabs( d[k*nX+i] - row( i )[k] ) * pval[i] / s );
    }
    return w; };
  auto zero_except = [&]( std::vector<double> const& d, std::vector<size_t> const& keep ){
    double z = 0.;
    for( size_t k = 0; k < nf; ++k ) for( size_t c = 0; c < nX; ++c )
      if( std::find( keep.begin(), keep.end(), c ) == keep.end() ) z = std::max( z, std::fabs( d[k*nX+c] ) );
    return z; };
  try{
    auto const dN = jac( {} );
    { std::ostringstream o; double const w = col_err( dN, { R_f, R_l, R_th } );
      o << "B2 numeric dF/d(R_f, R_l, R_th) == forward sensitivity (worst scaled " << w << ")"; check( w < 1e-5, o.str() ); }   // two FSA runs, 3 vs 10 directions: integrator tolerance
    { std::ostringstream o; double const z = zero_except( dN, { R_f, R_l, R_th } );
      o << "B2 numeric: 7 map-2 parameters and 5 constants ZERO by contract (max |.| " << z << ")"; check( z == 0., o.str() ); }
    auto const dS = jac( { cv[K], pv[L_d] } );
    { std::ostringstream o; double const w = col_err( dS, { L_d } );
      o << "B3 SYMDIFF dF/dL_d (map-2 parameter) == forward sensitivity (worst scaled " << w << ")"; check( w < 1e-6, o.str() ); }
    { // FD instrument: named solves of the declared model at RTOL = ATOL = 1e-13 (F ~ 1e-4, so the model's
      // ATOL 1e-10 would leave ~1e-10/(2h) of noise -- 1e-4 relative); h = 1e-4 K.  The op keeps its COPY.
      double const h = 1e-4 * cval[K];
      auto const oA = I.options.ATOL, oR = I.options.RTOL;  I.options.ATOL = I.options.RTOL = 1e-13;
      auto vIp = vIn, vIm = vIn;  vIp[np+K].vals[0] += h;  vIm[np+K].vals[0] -= h;
      I.solve( vIp );  auto const Fp = I.val_function();  I.solve( vIm );  auto const Fm = I.val_function();
      I.options.ATOL = oA;  I.options.RTOL = oR;
      double w = 0.;
      for( size_t k = 0; k < nf; ++k ){ double const fd = ( Fp[k] - Fm[k] ) / ( 2.*h );
        w = std::max( w, std::fabs( dS[k*nX+np+K] - fd ) / std::max( std::fabs( fd ), 1e-300 ) ); }
      std::ostringstream o; o << "B3 SYMDIFF dF/dK (CONSTANT, control gain) == central FD, solves at 1e-13 (worst rel " << w << ")"; check( w < 1e-5, o.str() ); }
    { std::ostringstream o; double const z = zero_except( dS, { L_d, np+K } );
      o << "B3 SYMDIFF: every other column ZERO by contract (max |.| " << z << ")"; check( z == 0., o.str() ); }
    auto const d1 = jac( { pv[R_f], pv[R_l], pv[R_th] } );              // the SECOND SYMDIFF on this op
    { std::ostringstream o; double const w = col_err( d1, { R_f, R_l, R_th } );
      o << "B4 SYMDIFF on map 1 (second SYMDIFF on the same op) == forward sensitivity (worst scaled " << w << ")";
      check( w < 1e-5, o.str() ); }                                                // augmented system vs reference run: integrator tolerance
  } catch( std::exception& e ){ check( false, std::string( "B2-B4 threw: " ) + e.what() ); }

  std::cout << "\n  ODESLV_NIFTE: " << npass << " passed, " << nfail << " failed -- " << ( nfail? "FAILURES": "ALL PASS" ) << "\n";
  return nfail? 1: 0;
}
