// DAESLV_moltime -- MoL phase 6c, time-varying inputs: a model distributed in space with inputs that vary in time, on
// DAESLV (each spatial slot an input over t of the semi-discrete DAE, with the declared collocation in t; the levels
// the integrator's parameters) and on OCFESLV, declared once.
// (A) a time-varying INLET: u_t = alpha u_zz, u(0,t) = g(t) with g an input over {t}, 3 LGL nodes per stage on 2
//     stages, its levels sampled from a smooth profile so g is continuous (OCFESLV collocates t with the state
//     continuous at the stage boundary, so a jump in Dirichlet data is not a case it can take), u_z(1,t) = 0,
//     u0 = g(t0); outputs the inventory at T and the outlet integral.  Values against OCFESLV (a ladder in its time
//     order); d/dg_k forward = adjoint = central difference, including through the INITIAL row u0 = g(t0) -- which
//     ic_sens lost until 5 Nov 2026 (the placeholder chain; also fixed for the ODESLV forward path).
// (B) a DISTRIBUTED, time-varying SOURCE s(t,z) with its own collocation, 3 LGL nodes per element in t AND z (2 stages,
//     2 spatial elements): 36 declared DOFs whose input_dofs order interleaves t and z -- the bijection declared DOF
//     <-> working parameter, since the semi-discretisation keeps one input over t per spatial slot.  J = int u dz at
//     T against OCFESLV (the same spatial discretisation: a source at a boundary or interface node enters through the
//     row structure there, so no quadrature-weight closed form holds); forward = adjoint = central difference, per DOF.
// (D) a TRANSITION on a field (MoL 6c): the zero-flux rod with a source p e^-t, u0 = 1, and u(tau+) = u(tau-)/2 at
//     tau = 0.5, an element boundary of t.  Applied at the field's differential nodes, both sides composed there; the
//     integrator's stop machinery re-solves the algebraic nodes and carries the sensitivities across.  A uniform field
//     is in the operator's kernel, so with J = int u dz at T: J = [1 + p(1 - e^-tau)]/2 + p(e^-tau - e^-T) and
//     dJ/dp = (1 - e^-tau)/2 + e^-tau - e^-T, exactly on any mesh -- value, forward and adjoint against them.
//     (D2) the same at tau = 0.3, INSIDE element 0 of t: a stop of its own, as ODESLV (refused until 12 Nov 2026).
#include <iostream>
#include <iomanip>
#include <cmath>
#include "ocfeslv.hpp"
#include "daeslv_idas.hpp"
#include "odeslv_cvodes.hpp"
using namespace mc;

static int nfail = 0;
static void check( std::string const& w, bool ok )
{ std::cout << "  " << std::left << std::setw(88) << w << ( ok? "PASS": "FAIL" ) << "\n"; if( !ok ) ++nfail; }
static double rel( double a, double b ){ return std::fabs( a - b )/std::max( 1e-300, std::fabs( b ) ); }
static double const ALPHA = 0.1, TEND = 1.;

static void colloc_options( OCFESLV& oc )
{ oc.options.INTERFACE.IMPOSITION = OCFESLV::Options::IC_STRONG;  oc.options.DISPLAY_LEVEL = 0;  oc.options.SOLVE.MARCHING = false; }
static void colloc_options( DAESLV_IDAS& S )
{ S.options.INTERFACE.IMPOSITION = OCFEMOL::Options::IC_STRONG;  S.options.DISPLAY_LEVEL = 0;
  S.options.SOLVE.RTOL = 1e-10;  S.options.SOLVE.ATOL = 1e-12;  S.options.SOLVE.RTOLB = 1e-10;  S.options.SOLVE.ATOLB = 1e-12;
  S.options.SOLVE.VERBOSE = 0; }

// (A) the inlet rod; NODT: the time order (OCFESLV collocates t; DAESLV takes the t domain's elements as its stages)
template <typename OC>
static void declare_inlet( OC& oc, FFGraph& G, FFVar& t, FFVar& z, FFVar& u, FFVar& g, size_t NODT )
{
  t = G.add_var("t"); z = G.add_var("z"); u = G.add_var("u"); g = G.add_var("g");
  oc.add_domain( t, FFDom( 0., TEND, 2, FFDom::LGR, NODT ) );
  oc.add_domain( z, FFDom( 0., 1., 3, FFDom::LGL, 5 ) );
  oc.set_evolution_domain( t );
  oc.add_state( u, {t,z} ); oc.update_ref( u, 1. );
  oc.add_input( g, {t}, FFDom::LGL, 3 );                              // 3 levels per stage, continuous across them
  FFPartial OpP;  FFIntegral OpI;
  int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
  typename OC::EqnOptions io( OC::EqnRole::INTERIOR, 0 ), ib( OC::EqnRole::BOUNDARY, 0 ), ii( OC::EqnRole::INITIAL, 0 );
  oc.add_equation( OpP(u,t) - OpP( ALPHA*OpP(u,z), z ), {t,z}, {TI, ZI},                io );
  oc.add_equation( u - g,                              {t,z}, {TI, FFDom::LB},         ib );
  oc.add_equation( OpP(u,z),                           {t,z}, {TI, FFDom::UB},         ib );
  oc.add_equation( u - g,                              {t,z}, {FFDom::LB, FFDom::ALL}, ii );
  oc.add_output( OpI( u, z ), {t}, {TEND} );                         // 0: the inventory at T
  oc.add_output( OpI( u, t ), {z}, {1.} );                           // 1: the outlet integral
  colloc_options( oc );
}

static bool inlet_daeslv( int mode, std::vector<double> const& G, std::vector<double>& f, std::vector<std::vector<double>>& Gr, std::string& why )
{
  FFGraph dag;  DAESLV_IDAS S( &dag );  FFVar t, z, u, g;
  declare_inlet( S, dag, t, z, u, g, 3 );
  if( !S.setup() ){ why = "setup: " + S.extract_error(); return false; }
  std::vector<FFModel::InputVal> const in{ { g, G } };
  auto const st = mode == 0? S.solve( in ): mode == 1? S.solve_fsens( in ): S.solve_asens( in );
  if( st != DAESLV_IDAS::NORMAL ){ why = "solve: " + S.extract_error(); return false; }
  f = S.val_function();  if( mode ) Gr = S.val_function_gradient();
  return f.size() == 2;
}
static bool inlet_ocfeslv( std::vector<double> const& G, size_t NODT, std::vector<double>& f, std::string& why )
{
  FFGraph dag;  OCFESLV oc( &dag );  FFVar t, z, u, g;
  declare_inlet( oc, dag, t, z, u, g, NODT );
  if( !oc.setup() ){ why = "setup"; return false; }
  if( oc.solve( std::vector<FFModel::InputVal>{ { g, G } } ) != OCFESLV::NORMAL ){ why = "solve"; return false; }
  f = oc.val_function();  return f.size() == 2;
}

// (B) the source rod
template <typename OC>
static void declare_source( OC& oc, FFGraph& G, FFVar& t, FFVar& z, FFVar& u, FFVar& s, size_t NODT )
{
  t = G.add_var("t"); z = G.add_var("z"); u = G.add_var("u"); s = G.add_var("s");
  oc.add_domain( t, FFDom( 0., TEND, 2, FFDom::LGR, NODT ) );
  oc.add_domain( z, FFDom( 0., 1., 2, FFDom::LGL, 5 ) );
  oc.set_evolution_domain( t );
  oc.add_state( u, {t,z} ); oc.update_ref( u, 1. );
  oc.add_input( s, {t,z}, FFDom::LGL, 3 );                            // 3 LGL nodes per element, in t and in z
  FFPartial OpP;  FFIntegral OpI;
  int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
  typename OC::EqnOptions io( OC::EqnRole::INTERIOR, 0 ), ib( OC::EqnRole::BOUNDARY, 0 ), ii( OC::EqnRole::INITIAL, 0 );
  oc.add_equation( OpP(u,t) - OpP( ALPHA*OpP(u,z), z ) - s, {t,z}, {TI, ZI},                io );
  oc.add_equation( OpP(u,z),                               {t,z}, {TI, FFDom::LB},         ib );
  oc.add_equation( OpP(u,z),                               {t,z}, {TI, FFDom::UB},         ib );
  oc.add_equation( u - 1.,                                 {t,z}, {FFDom::LB, FFDom::ALL}, ii );
  oc.add_output( OpI( u, z ), {t}, {TEND} );
  colloc_options( oc );
}

int main()
{
  std::set_terminate( []{ try{ throw; } catch( mc::FFBase::Exceptions& e ){ std::cerr << "FFBase::Exceptions: " << e.what() << "\n"; }
                          catch( std::exception const& e ){ std::cerr << "exception: " << e.what() << "\n"; } catch( ... ){} std::abort(); } );
  std::cout << std::scientific << std::setprecision(3);
  std::string why;
  // ---- (A)
  std::cout << "(A) a time-varying inlet g(t), piecewise constant on 2 stages\n";
  // g at the 3 LGL nodes of each stage ([0, .5], [.5, 1]): 1 + 0.5 sin(pi t) -- the shared node gets one value
  double const tn[6] = { 0., 0.25, 0.5, 0.5, 0.75, 1. };
  std::vector<double> G( 6 );  for( int k = 0; k < 6; ++k ) G[k] = 1. + 0.5*std::sin( M_PI*tn[k] );
  double const H = 1e-2;                                                // the model is linear in g
  std::vector<double> f0, fF, fA, fp, fm;  std::vector<std::vector<double>> gF, gA, d;
  bool const v = inlet_daeslv( 0, G, f0, d, why );
  check( "(A) DAESLV: the model with a time-varying input sets up and solves by name" + ( v? std::string(): " -- " + why ), v );
  bool const fw = inlet_daeslv( 1, G, fF, gF, why );
  check( "(A) forward sensitivities" + ( fw? std::string(): " -- " + why ), fw && gF.size() == 6 && gF[0].size() == 2 );
  bool const ad = inlet_daeslv( 2, G, fA, gA, why );
  check( "(A) adjoint sensitivities" + ( ad? std::string(): " -- " + why ), ad && gA.size() == 6 && gA[0].size() == 2 );
  if( v && fw && ad ){
    // errors scaled by each output's largest gradient (a near-zero entry has no relative error of its own)
    double gmax[2] = { 0., 0. };  for( size_t k = 0; k < 6; ++k ) for( int i = 0; i < 2; ++i ) gmax[i] = std::max( gmax[i], std::fabs( gF[k][i] ) );
    double eFA = 0., eFD = 0.;
    for( size_t k = 0; k < 6; ++k ){
      std::vector<double> Gp = G, Gm = G;  Gp[k] += H;  Gm[k] -= H;
      inlet_daeslv( 0, Gp, fp, d, why );  inlet_daeslv( 0, Gm, fm, d, why );
      for( int i = 0; i < 2; ++i ){
        double const fd = ( fp[i] - fm[i] )/( 2.*H );               // the model is linear in g: no truncation error
        eFA = std::max( eFA, std::fabs( gA[k][i] - gF[k][i] )/gmax[i] );  eFD = std::max( eFD, std::fabs( gF[k][i] - fd )/gmax[i] );
        if( k < 2 || i == 0 ) std::cout << "    d f" << i << " / d g" << k << ":  fwd " << gF[k][i] << "  adj " << gA[k][i] << "  FD " << fd << "\n";
      }
    }
    std::cout << "    max scaled: adj-fwd " << eFA << "  fwd-FD " << eFD << "\n";
    check( "(A) adjoint = forward, to 1e-8 of the largest gradient",             eFA <= 1e-8 );
    check( "(A) forward = central difference, to 1e-6 of the largest gradient",  eFD <= 1e-6 );
    std::vector<std::vector<double>> fL( 3 );  bool oc = true;
    for( int k = 0; k < 3 && oc; ++k ) oc = inlet_ocfeslv( G, 3 + k, fL[k], why );
    check( "(A) OCFESLV at time orders 3, 4, 5" + ( oc? std::string(): " -- " + why ), oc );
    if( oc ){
      bool conv = true;
      for( int i = 0; i < 2; ++i ){
        std::cout << "    f" << i << "  DAESLV " << std::setprecision(10) << f0[i] << std::setprecision(2) << "   |DAESLV - OCFESLV|/|.| orders 3 4 5:";
        for( int k = 0; k < 3; ++k ) std::cout << " " << rel( f0[i], fL[k][i] );
        std::cout << std::setprecision(3) << "\n";
        for( int k = 0; k < 2; ++k ) conv = conv && rel( f0[i], fL[k][i] ) >= 3.*rel( f0[i], fL[k+1][i] );
      }
      check( "(A) the gap to OCFESLV falls >= 3x per time order", conv );
      check( "(A) at order 5: DAESLV = OCFESLV to 1e-5",          rel( f0[0], fL[2][0] ) <= 1e-5 && rel( f0[1], fL[2][1] ) <= 1e-5 );
    }
  }
  // ---- (B)
  std::cout << "(B) a distributed, time-varying source s(t,z): (2 stages x 3 nodes) x (2 elements x 3 nodes) = 36 DOFs\n";
  {
    FFGraph dag;  DAESLV_IDAS S( &dag );  FFVar t, z, u, s;
    declare_source( S, dag, t, z, u, s, 3 );
    bool ok = S.setup();
    check( "(B) setup" + ( ok? std::string(): " -- " + S.extract_error() ), ok );
    if( ok ){
      auto const dofs = S.input_dofs( s );
      std::cout << "    declared DOFs: " << dofs.size() << "\n";
      check( "(B) 36 declared DOFs", dofs.size() == 36 );
      std::vector<double> SV( dofs.size() );                          // a distinct value per DOF
      for( size_t k = 0; k < dofs.size(); ++k ) SV[k] = 0.3 + 0.02*k;
      std::vector<FFModel::InputVal> const in{ { s, SV } };
      ok = S.solve( in ) == DAESLV_IDAS::NORMAL;
      check( "(B) value solve by name" + ( ok? std::string(): " -- " + S.extract_error() ), ok );
      if( ok ){
        double const J0 = S.val_function().at( 0 );
        std::vector<std::vector<double>> gF, gA;
        bool const fw = S.solve_fsens( in ) == DAESLV_IDAS::NORMAL;  if( fw ) gF = S.val_function_gradient();
        bool const ad = S.solve_asens( in ) == DAESLV_IDAS::NORMAL;  if( ad ) gA = S.val_function_gradient();
        check( "(B) forward and adjoint sensitivities (36 rows)", fw && ad && gF.size() == 36 && gA.size() == 36 );
        if( fw && ad && gF.size() == 36 && gA.size() == 36 ){
          double eFA = 0., eFD = 0.;
          for( size_t k = 0; k < 36; ++k ){
            auto const itn = dofs[k].node.find( z ), itt = dofs[k].node.find( t );
            std::vector<double> Sp = SV, Sm = SV;  Sp[k] += 0.1;  Sm[k] -= 0.1;   // J is linear in s: a large step
            S.solve( std::vector<FFModel::InputVal>{ { s, Sp } } );  double const Jp = S.val_function().at( 0 );
            S.solve( std::vector<FFModel::InputVal>{ { s, Sm } } );  double const Jm = S.val_function().at( 0 );
            double const fd = ( Jp - Jm )/0.2;
            eFA = std::max( eFA, rel( gA[k][0], gF[k][0] ) );  eFD = std::max( eFD, rel( gF[k][0], fd ) );
            if( k < 4 ) std::cout << "    DOF " << k << " (t el " << dofs[k].element.at( t ) << " node " << itt->second << ", z el " << dofs[k].element.at( z ) << " node " << itn->second
                                  << "): fwd " << gF[k][0] << " adj " << gA[k][0] << " FD " << fd << "\n";
          }
          std::cout << "    J " << J0 << "   max rel: adj-fwd " << eFA << "  fwd-FD " << eFD << "\n";
          check( "(B) adjoint = forward, to 1e-8",                         eFA <= 1e-8 );
          check( "(B) forward = central difference, to 1e-6",              eFD <= 1e-6 );
          // the value against OCFESLV, on the same declaration (orders 3, 4, 5)
          std::vector<double> fo( 3, 0. );  bool oc = true;
          for( int k = 0; k < 3 && oc; ++k ){
            FFGraph dag2;  OCFESLV oc2( &dag2 );  FFVar t2, z2, u2, s2;
            declare_source( oc2, dag2, t2, z2, u2, s2, 3 + k );
            oc = oc2.setup() && oc2.solve( std::vector<FFModel::InputVal>{ { s2, SV } } ) == OCFESLV::NORMAL;
            if( oc ) fo[k] = oc2.val_function().at( 0 );
          }
          std::cout << "    |J - OCFESLV|/|.| orders 3 4 5: " << rel( J0, fo[0] ) << " " << rel( J0, fo[1] ) << " " << rel( J0, fo[2] ) << "\n";
          check( "(B) OCFESLV solves the same declaration",                 oc );
          check( "(B) at order 5: J = OCFESLV to 1e-5",                     oc && rel( J0, fo[2] ) <= 1e-5 );
        }
      }
    }
  }
  // ---- (C) the gap that (A) exposed, in its simplest form, on both integrators
  std::cout << "(C) an INITIAL row reading a time-varying input: x' = -x, x(0) = u(t0), u piecewise constant on 2 stages\n";
  {
    double const e1 = std::exp( -1. );
    auto ode = [&]( auto& s, FFGraph& G, int mode, double& x1, double& du0, double& du1, std::string& w ){
      FFVar t = G.add_var("t"), x = G.add_var("x"), u = G.add_var("u");
      s.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 3 ) );  s.set_evolution_domain( t );
      s.add_state( x, {t} );  s.update_ref( x, 1. );
      s.add_input( u, {t}, FFDom::LGL, 1 );
      FFPartial OpP;  FFEval OpE;
      using OC = std::decay_t<decltype( s )>;
      typename OC::EqnOptions io( OC::EqnRole::INTERIOR, 0 ), ii( OC::EqnRole::INITIAL, 0 );
      s.add_equation( OpP(x,t) + x, {t}, {FFDom::ALL - FFDom::LB}, io );
      s.add_equation( x - u,        {t}, {FFDom::LB},              ii );
      s.add_output( OpE( x, t, 1. ) );
      s.options.SOLVE.RTOL = 1e-10;  s.options.SOLVE.ATOL = 1e-12;  s.options.SOLVE.RTOLB = 1e-10;  s.options.SOLVE.ATOLB = 1e-12;
      s.options.DISPLAY_LEVEL = 0;  s.options.SOLVE.VERBOSE = 0;
      if( !s.setup() ){ w = "setup"; return false; }
      std::vector<FFModel::InputVal> const in{ { u, { 0.8, 1.5 } } };
      if( ( mode == 1? s.solve_fsens( in ): s.solve_asens( in ) ) != OC::NORMAL ){ w = "solve: " + s.extract_error(); return false; }
      x1 = s.val_function().at(0);  du0 = s.val_function_gradient().at(0).at(0);  du1 = s.val_function_gradient().at(1).at(0);
      return true; };
    for( int mode : { 1, 2 } ){
      double x1 = 0., a = 0., b = 0.;  std::string w;
      { FFGraph G;  DAESLV_IDAS S( &G );  bool const ok = ode( S, G, mode, x1, a, b, w );
        check( std::string( "(C) DAESLV " ) + ( mode == 1? "forward": "adjoint" ) + ": dx(1)/du0 = e^-1, dx(1)/du1 = 0" + ( ok? "": " -- " + w ),
               ok && rel( x1, 0.8*e1 ) <= 1e-8 && rel( a, e1 ) <= 1e-8 && std::fabs( b ) <= 1e-12 ); }
      { FFGraph G;  ODESLV_CVODES S( &G );  bool const ok = ode( S, G, mode, x1, a, b, w );
        bool const good = ok && rel( x1, 0.8*e1 ) <= 1e-8 && rel( a, e1 ) <= 1e-8 && std::fabs( b ) <= 1e-12;
        if( mode == 1 ) check( "(C) ODESLV forward: the same (the explicit bypass yields to the solve)" + ( ok? std::string(): " -- " + w ), good );
        else check( "(C) ODESLV adjoint: the same (its initial-condition chain through ic_sens, as DAESLV's)" + ( ok? std::string(): " -- " + w ), good ); }
    }
  }
  // ---- (H) a transition reading an input DISTRIBUTED IN SPACE (10 Nov 2026; refused before): the zero-flux rod with
  //      u0 = 1 and source p e^-t, u(tau+) = u(tau-) + w(z), w on its OWN coarser collocation (2 LGL nodes per element
  //      against the field's 5), so it is read at the field's nodes through its Lagrange interpolant.  (H1) a POINTWISE
  //      field, u_t = p e^-t (1 + z) (z read, as FFModel requires of an equation declared over it): u (T, z) is linear
  //      per element, the field's quadrature exact for it -- int u dz (T) = 1 + 3/2 p (1 - e^-T) + int w dz EXACTLY on
  //      the mesh, int w dz = sum_e h/2 (w_l + w_r), d/dw_k = h/2 for each of w's 6 DOFs.  (H2) the diffusive rod: strong-form collocation conserves the discrete integral only to ~2e-4 once
  //      the field is not uniform, so its gate is forward = adjoint.
  std::cout << "(H) a transition reading an input distributed in space: u(tau+) = u(tau-) + w(z), w on 2 LGL nodes per element\n";
  {
    double const P = 0.7, TE = 1., TAU = 0.5;
    std::vector<double> const wv{ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 };
    double const h = 1./3., Iw = 0.5*h*( wv[0]+wv[1] + wv[2]+wv[3] + wv[4]+wv[5] );
    for( int diffusive = 0; diffusive < 2; ++diffusive ){
    double const Je = 1. + ( diffusive? 1.: 1.5 )*P*( 1. - std::exp( -TE ) ) + Iw;
    std::vector<double> F[3];  std::vector< std::vector<double> > Gr[3];  bool ok = true;  std::string w;
    std::string const tag = diffusive? "(H2) the diffusive rod": "(H1) a pointwise field";
    for( int mode = 0; mode < 3 && ok; ++mode ){
      FFGraph G;  DAESLV_IDAS S( &G );
      FFVar t = G.add_var( "t" ), z = G.add_var( "z" ), u = G.add_var( "u" ), p = G.add_var( "p" ), wz = G.add_var( "w" );
      S.add_domain( t, FFDom( 0., TE, 2, FFDom::LGR, 3 ) );  S.add_domain( z, FFDom( 0., 1., 3, FFDom::LGL, 5 ) );
      S.set_evolution_domain( t );  S.add_state( u, {t,z} );  S.update_ref( u, 1. );  S.add_input( p, {} );
      S.add_input( wz, {z}, FFDom::LGL, 2 );
      FFPartial OpP;  FFIntegral OpI;
      int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
      DAESLV_IDAS::EqnOptions io( DAESLV_IDAS::EqnRole::INTERIOR, 0 ), ib( DAESLV_IDAS::EqnRole::BOUNDARY, 0 ), ii( DAESLV_IDAS::EqnRole::INITIAL, 0 );
      if( diffusive ){
        S.add_equation( OpP(u,t) - OpP( 0.1*OpP(u,z), z ) - p*exp( -t ), {t,z}, {TI, ZI}, io );
        S.add_equation( OpP(u,z), {t,z}, {TI, FFDom::LB}, ib );  S.add_equation( OpP(u,z), {t,z}, {TI, FFDom::UB}, ib );
      }
      else S.add_equation( OpP(u,t) - p*exp( -t )*( 1. + z ), {t,z}, {TI, FFDom::ALL}, io );
      S.add_equation( u - 1., {t,z}, {FFDom::LB, FFDom::ALL}, ii );
      S.add_transition( { u + wz }, { u }, t, TAU );
      S.add_output( OpI( u, z ), {t}, {TE} );
      S.options.SOLVE.RTOL = S.options.SOLVE.ATOL = 1e-11;  S.options.SOLVE.VERBOSE = 0;  S.options.DISPLAY_LEVEL = 0;
      if( !S.setup() ){ ok = false;  w = "setup: " + S.extract_error();  break; }
      std::vector<FFModel::InputVal> in{ { p, { P } }, { wz, wv } };
      auto const st = mode == 0? S.solve( in ): mode == 1? S.solve_fsens( in ): S.solve_asens( in );
      if( st != DAESLV_IDAS::NORMAL ){ ok = false;  w = "solve: " + S.extract_error();  break; }
      F[mode] = S.val_function();  if( mode ) Gr[mode] = S.val_function_gradient();
    }
    check( tag + ": set up and solved, value, forward and adjoint (refused before)" + ( ok? std::string(): " -- " + w ), ok );
    double ev = NAN, ef = 0., ea = 0., efa = 0.;  size_t nw = 0;
    if( ok ){
      ev = rel( F[0].at(0), Je );
      // the gradient rows: p first, then w's DOFs in declaration order
      for( size_t k = 1; k < Gr[1].size() && k < Gr[2].size(); ++k ){
        ef = std::max( ef, rel( Gr[1][k].at(0), 0.5*h ) );  ea = std::max( ea, rel( Gr[2][k].at(0), 0.5*h ) );
        efa = std::max( efa, rel( Gr[2][k].at(0), Gr[1][k].at(0) ) );  ++nw; }
      std::cout << std::scientific << std::setprecision(3) << "    " << tag << ": int u dz (T) rel err " << ev << ";  " << nw
                << " w DOFs, d/dw_k vs h/2: forward " << ef << ", adjoint " << ea << ";  adjoint vs forward " << efa << "\n";
    }
    if( !diffusive )
      check( tag + ": int u dz (T) and d/dw_k = h/2 for all 6 DOFs = closed forms, value and forward to 1e-8, adjoint to 1e-6",
             ok && nw == 6 && ev < 1e-8 && ef < 1e-8 && ea < 1e-6 );
    else
      check( tag + ": the adjoint gradient = the forward one to 1e-6 for all 6 DOFs (value within 1e-3: collocation's conservation)",
             ok && nw == 6 && efa < 1e-6 && ev < 1e-3 );
    }
  }

  // ---- (D), and (D2) off the element boundaries of t
  for( double const TAU : { 0.5, 0.3 } ){
  std::string const D = TAU == 0.5? "(D)": "(D2)";
  if( TAU == 0.5 ) std::cout << "(D) a transition on a field: u(tau+) = u(tau-)/2 at tau = 0.5\n";
  else             std::cout << "(D2) the same transition at tau = 0.3, inside element 0 of t (a stop of its own, as ODESLV)\n";
  {
    double const P = 0.7, eT = std::exp( -TEND ), eTau = std::exp( -TAU );
    double const Jex = 0.5*( 1. + P*( 1. - eTau ) ) + P*( eTau - eT ), dJex = 0.5*( 1. - eTau ) + eTau - eT;
    auto rod = [&]( int mode, double& J, double& dJ, std::string& w ){
      FFGraph G;  DAESLV_IDAS S( &G );
      FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u"), p = G.add_var("p");
      S.add_domain( t, FFDom( 0., TEND, 2, FFDom::LGR, 3 ) );
      S.add_domain( z, FFDom( 0., 1., 3, FFDom::LGL, 5 ) );
      S.set_evolution_domain( t );
      S.add_state( u, {t,z} ); S.update_ref( u, 1. );
      S.add_input( p, {} );
      FFPartial OpP;  FFIntegral OpI;
      int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
      DAESLV_IDAS::EqnOptions io( DAESLV_IDAS::EqnRole::INTERIOR, 0 ), ib( DAESLV_IDAS::EqnRole::BOUNDARY, 0 ), ii( DAESLV_IDAS::EqnRole::INITIAL, 0 );
      S.add_equation( OpP(u,t) - OpP( ALPHA*OpP(u,z), z ) - p*exp( -t ), {t,z}, {TI, ZI},                io );
      S.add_equation( OpP(u,z),                                        {t,z}, {TI, FFDom::LB},         ib );
      S.add_equation( OpP(u,z),                                        {t,z}, {TI, FFDom::UB},         ib );
      S.add_equation( u - 1.,                                          {t,z}, {FFDom::LB, FFDom::ALL}, ii );
      S.add_transition( { 0.5*u }, { u }, t, TAU );
      S.add_output( OpI( u, z ), {t}, {TEND} );
      colloc_options( S );
      if( !S.setup() ){ w = "setup: " + S.extract_error(); return false; }
      std::vector<FFModel::InputVal> const in{ { p, { P } } };
      auto const st = mode == 0? S.solve( in ): mode == 1? S.solve_fsens( in ): S.solve_asens( in );
      if( st != DAESLV_IDAS::NORMAL ){ w = "solve: " + S.extract_error(); return false; }
      J = S.val_function().at( 0 );  if( mode ) dJ = S.val_function_gradient().at( 0 ).at( 0 );
      return true; };
    double J0 = 0., Jf = 0., Ja = 0., dF = 0., dA = 0., d0 = 0.;  std::string w;
    bool const rv = rod( 0, J0, d0, w );
    check( D + " a model with a transition sets up and solves" + ( rv? std::string(): " -- " + w ), rv );
    bool const rf = rod( 1, Jf, dF, w );
    check( D + " forward sensitivities across the stop" + ( rf? std::string(): " -- " + w ), rf );
    bool const ra = rod( 2, Ja, dA, w );
    check( D + " adjoint sensitivities across the stop" + ( ra? std::string(): " -- " + w ), ra );
    if( rv && rf && ra ){
      std::cout << "    J " << J0 << " (closed form " << Jex << ")   dJ/dp: fwd " << dF << " adj " << dA << " (closed form " << dJex << ")\n";
      check( D + " J = closed form, to 1e-8",                          rel( J0, Jex ) <= 1e-8 );
      check( D + " dJ/dp = closed form, forward and adjoint, to 1e-8", rel( dF, dJex ) <= 1e-8 && rel( dA, dJex ) <= 1e-8 );
    }
  }
  }
  // ---- (E) NONLINEAR functions of horizon integrals (7 Nov 2026).  The zero-flux rod with u0 = 1 and the source
  //      p e^-t: a uniform field is in the operator's kernel, so J(t) = int u dz = 1 + p (1 - e^-t) exactly on any
  //      mesh, and A = int J dt, B = int J^2 dt, sqrt A, A/B have closed forms.  B's integrand is the SQUARE of a
  //      spatial integral: a reduction nested inside a nonlinearity, which the collocation's root-level assembly
  //      sums per element (J^2/3 on three elements, on OCFESLV too -- an open item there); under MoL that reduction
  //      is a scalar part of its own and the square is composed over it at every t (OCFEMOL::t_OutputTime).
  std::cout << "(E) nonlinear functions of horizon integrals: A, sqrt A, A/B with B = int (int u dz)^2 dt\n";
  {
    double const P = 0.7, e1 = std::exp( -TEND ), e2 = std::exp( -2.*TEND );
    double const Ae = 1. + P*e1, dAe = e1;
    double const Be = (1+P)*(1+P) - 2*P*(1+P)*(1-e1) + P*P*(1-e2)/2., dBe = 2*(1+P) - 2*(1+2*P)*(1-e1) + P*(1-e2);
    double const Fe[3] = { Ae, std::sqrt( Ae ), Ae/Be }, dFe[3] = { dAe, dAe/( 2*std::sqrt( Ae ) ), ( dAe*Be - Ae*dBe )/( Be*Be ) };
    auto rod = [&]( int mode, double* F, double* dF, std::string& w ){
      FFGraph G;  DAESLV_IDAS S( &G );
      FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u"), p = G.add_var("p");
      S.add_domain( t, FFDom( 0., TEND, 2, FFDom::LGR, 3 ) );
      S.add_domain( z, FFDom( 0., 1., 3, FFDom::LGL, 5 ) );
      S.set_evolution_domain( t );
      S.add_state( u, {t,z} ); S.update_ref( u, 1. );
      S.add_input( p, {} );
      FFPartial OpP;  FFIntegral OpI;
      int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
      DAESLV_IDAS::EqnOptions io( DAESLV_IDAS::EqnRole::INTERIOR, 0 ), ib( DAESLV_IDAS::EqnRole::BOUNDARY, 0 ), ii( DAESLV_IDAS::EqnRole::INITIAL, 0 );
      S.add_equation( OpP(u,t) - OpP( ALPHA*OpP(u,z), z ) - p*exp( -t ), {t,z}, {TI, ZI},                io );
      S.add_equation( OpP(u,z),                                        {t,z}, {TI, FFDom::LB},         ib );
      S.add_equation( OpP(u,z),                                        {t,z}, {TI, FFDom::UB},         ib );
      S.add_equation( u - 1.,                                          {t,z}, {FFDom::LB, FFDom::ALL}, ii );
      FFVar J = OpI( u, z ), A = OpI( J, t ), B = OpI( J*J, t );
      S.add_output( A );  S.add_output( sqrt( A ) );  S.add_output( A / B );
      colloc_options( S );
      if( !S.setup() ){ w = "setup: " + S.extract_error(); return false; }
      std::vector<FFModel::InputVal> const in{ { p, { P } } };
      auto const st = mode == 0? S.solve( in ): mode == 1? S.solve_fsens( in ): S.solve_asens( in );
      if( st != DAESLV_IDAS::NORMAL ){ w = "solve: " + S.extract_error(); return false; }
      for( int i = 0; i < 3; ++i ){ F[i] = S.val_function().at( i );  if( mode ) dF[i] = S.val_function_gradient().at( 0 ).at( i ); }
      return S.nf() == 3; };
    double F0[3], Ff[3], Fa[3], dF[3], dA[3], d0[3];  std::string w;
    bool const rv = rod( 0, F0, d0, w ), rf = rod( 1, Ff, dF, w ), ra = rod( 2, Fa, dA, w );
    check( "(E) three nonlinear outputs set up and solve, value, forward and adjoint" + ( rv && rf && ra? std::string(): " -- " + w ), rv && rf && ra );
    if( rv && rf && ra ){
      double ev = 0., ef = 0., ea = 0.;
      for( int i = 0; i < 3; ++i ){ ev = std::max( ev, rel( F0[i], Fe[i] ) );  ef = std::max( ef, rel( dF[i], dFe[i] ) );  ea = std::max( ea, rel( dA[i], dFe[i] ) ); }
      std::cout << "    A " << F0[0] << " sqrtA " << F0[1] << " A/B " << F0[2] << "   d/dp fwd " << dF[0] << " " << dF[1] << " " << dF[2]
                << "   max rel: value " << ev << " fwd " << ef << " adj " << ea << "\n";
      check( "(E) values = closed forms, to 1e-8",                      ev < 1e-8 );
      check( "(E) forward gradients = closed forms, to 1e-8",           ef < 1e-8 );
      check( "(E) adjoint gradients = closed forms, to 1e-6",           ea < 1e-6 );
    }
  }
  // ---- (F) TRAJECTORIES and (G) TRANSITIONS READING INPUTS (7 Nov 2026), on the zero-flux rod with u0 = 1 and source
  //      p e^-t: u(t,z) = 1 + p (1 - e^-t) at every node, on any mesh.  (G) also on the pointwise ODE x' = 0, where the
  //      transition x(tau+) = x(tau-) + g read a time-varying g with NO adjoint gradient (ic_adjoint lacked ic_sens's
  //      placeholder term; ODESLV was right) -- found here, fixed in BASE_DE.
  std::cout << "(F) trajectories and (G) transitions reading inputs, on the rod with a uniform field\n";
  {
    double const P = 0.7, R = 0.6, TAU = 0.5, e5 = std::exp( -TAU ), e1 = std::exp( -TEND );
    auto rodF = [&]( int mode, int variant, std::vector<double>& F, std::vector<std::vector<double>>& Gr, std::string& w ){
      FFGraph G;  DAESLV_IDAS S( &G );
      FFVar t = G.add_var("t"), z = G.add_var("z"), u = G.add_var("u"), p = G.add_var("p"), r = G.add_var("r"), g = G.add_var("g");
      S.add_domain( t, FFDom( 0., TEND, 2, FFDom::LGR, 3 ) );  S.add_domain( z, FFDom( 0., 1., 3, FFDom::LGL, 5 ) );
      S.set_evolution_domain( t );  S.add_state( u, {t,z} );  S.update_ref( u, 1. );  S.add_input( p, {} );
      if( variant == 1 ) S.add_input( r, {} );
      if( variant == 2 ) S.add_input( g, {t}, FFDom::LGL, 1 );                   // one level per stage
      FFPartial OpP;  FFEval OpE;  FFIntegral OpI;
      int const TI = FFDom::ALL - FFDom::LB, ZI = FFDom::ALL - FFDom::LB - FFDom::UB;
      DAESLV_IDAS::EqnOptions io( DAESLV_IDAS::EqnRole::INTERIOR, 0 ), ib( DAESLV_IDAS::EqnRole::BOUNDARY, 0 ), ii( DAESLV_IDAS::EqnRole::INITIAL, 0 );
      S.add_equation( OpP(u,t) - OpP( ALPHA*OpP(u,z), z ) - p*exp( -t ), {t,z}, {TI, ZI}, io );
      S.add_equation( OpP(u,z), {t,z}, {TI, FFDom::LB}, ib );  S.add_equation( OpP(u,z), {t,z}, {TI, FFDom::UB}, ib );
      S.add_equation( u - 1., {t,z}, {FFDom::LB, FFDom::ALL}, ii );
      if( variant == 0 ){ S.add_output( OpE( u, z, 1. ), {t}, std::vector<int>{ FFDom::ALL } );
                          S.add_output( u, {t,z}, std::vector<int>{ FFDom::ALL - FFDom::LB, FFDom::ALL } ); }
      if( variant == 1 ){ S.add_transition( { r*u }, { u }, t, TAU );  S.add_output( OpI( u, z ), {t}, {TEND} ); }
      if( variant == 2 ){ S.add_transition( { u + g }, { u }, t, TAU );  S.add_output( OpI( u, z ), {t}, {TEND} ); }
      colloc_options( S );
      if( !S.setup() ){ w = "setup: " + S.extract_error(); return false; }
      std::vector<FFModel::InputVal> in{ { p, { P } } };
      if( variant == 1 ) in.push_back( { r, { R } } );
      if( variant == 2 ) in.push_back( { g, { 0.3, 0.8 } } );
      auto const st = mode == 0? S.solve( in ): mode == 1? S.solve_fsens( in ): S.solve_asens( in );
      if( st != DAESLV_IDAS::NORMAL ){ w = "solve: " + S.extract_error(); return false; }
      F = S.val_function();  if( mode ) Gr = S.val_function_gradient();
      return true; };
    for( int v = 0; v < 3; ++v ){
      std::vector<double> F0, Ff, Fa;  std::vector<std::vector<double>> d0, Gf, Ga;  std::string w;
      bool const ok = rodF( 0, v, F0, d0, w ) && rodF( 1, v, Ff, Gf, w ) && rodF( 2, v, Fa, Ga, w );
      char const* name[3] = { "(F) trajectories: u(t,1) at the 3 stage times, the field at the last 2",
                              "(G) a transition reading a parameter, u(tau+) = r u(tau-)",
                              "(G) a transition reading a time-varying input, u(tau+) = u(tau-) + g" };
      check( std::string( name[v] ) + ": set up and solved, value, forward and adjoint" + ( ok? std::string(): " -- " + w ), ok );
      if( !ok ) continue;
      double ev = 0., ef = 0., ea = 0.;
      if( v == 0 ){                                   // every row is 1 + p (1 - e^-t_k), d/dp = 1 - e^-t_k
        double const tk[3] = { 0., TAU, TEND };
        ev = F0.size() == 3 + 2*15? 0.: 1.;
        for( size_t i = 0; i < F0.size(); ++i ){
          double const tt = i < 3? tk[i]: tk[ 1 + ( i - 3 ) % 2 ];
          ev = std::max( ev, rel( F0[i], 1. + P*( 1. - std::exp( -tt ) ) ) );
          ef = std::max( ef, std::fabs( Gf[0][i] - ( 1. - std::exp( -tt ) ) ) );  ea = std::max( ea, std::fabs( Ga[0][i] - ( 1. - std::exp( -tt ) ) ) );
        }
      }
      else if( v == 1 ){                              // J = r (1 + p (1 - e^-tau)) + p (e^-tau - e^-T)
        ev = rel( F0[0], R*( 1. + P*( 1. - e5 ) ) + P*( e5 - e1 ) );
        ef = std::max( rel( Gf[0][0], R*( 1. - e5 ) + e5 - e1 ), rel( Gf[1][0], 1. + P*( 1. - e5 ) ) );
        ea = std::max( rel( Ga[0][0], R*( 1. - e5 ) + e5 - e1 ), rel( Ga[1][0], 1. + P*( 1. - e5 ) ) );
      }
      else{                                           // the LEFT level: J = 1 + p (1 - e^-T) + 0.3, dJ/dg = (1, 0)
        ev = rel( F0[0], 1. + P*( 1. - e1 ) + 0.3 );
        ef = std::max( std::fabs( Gf[1][0] - 1. ), std::fabs( Gf[2][0] ) );
        ea = std::max( std::fabs( Ga[1][0] - 1. ), std::fabs( Ga[2][0] ) );
      }
      std::cout << "    max error: value " << ev << "  forward " << ef << "  adjoint " << ea << "\n";
      check( std::string( name[v] ) + ": = closed forms, values and forward to 1e-8, adjoint to 1e-6", ev < 1e-8 && ef < 1e-8 && ea < 1e-6 );
    }
    // the pointwise twin of the time-varying case, on both integrators
    auto ode = [&]( auto* tag, std::vector<std::vector<double>>& Gf, std::vector<std::vector<double>>& Ga, double& J, std::string& w ){
      using S_t = std::remove_pointer_t<decltype( tag )>;
      FFGraph G;  S_t s( &G );  FFVar t = G.add_var("t"), x = G.add_var("x"), g = G.add_var("g");
      s.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 3 ) );  s.set_evolution_domain( t );
      s.add_state( x, {t} );  s.update_ref( x, 1. );  s.add_input( g, {t}, FFDom::LGL, 1 );
      FFPartial OpP;  FFEval OpE;
      typename S_t::EqnOptions io( S_t::EqnRole::INTERIOR, 0 ), ii( S_t::EqnRole::INITIAL, 0 );
      s.add_equation( OpP(x,t), {t}, {FFDom::ALL - FFDom::LB}, io );  s.add_equation( x - 1., {t}, {FFDom::LB}, ii );
      s.add_transition( { x + g }, { x }, t, 0.5 );  s.add_output( OpE( x, t, 1. ) );
      s.options.SOLVE.RTOL = 1e-10;  s.options.SOLVE.ATOL = 1e-12;  s.options.SOLVE.VERBOSE = 0;  s.options.DISPLAY_LEVEL = 0;
      if( !s.setup() ){ w = "setup: " + s.extract_error(); return false; }
      std::vector<FFModel::InputVal> in{ { g, { 0.3, 0.8 } } };
      if( s.solve_fsens( in ) != S_t::NORMAL ){ w = "fsens"; return false; }  Gf = s.val_function_gradient();  J = s.val_function()[0];
      if( s.solve_asens( in ) != S_t::NORMAL ){ w = "asens"; return false; }  Ga = s.val_function_gradient();
      return true; };
    for( int which = 0; which < 2; ++which ){
      std::vector<std::vector<double>> Gf, Ga;  double J = 0.;  std::string w;
      bool const ok = which == 0? ode( (DAESLV_IDAS*)nullptr, Gf, Ga, J, w ): ode( (ODESLV_CVODES*)nullptr, Gf, Ga, J, w );
      double const err = ok? std::max( { std::fabs( J - 1.3 ), std::fabs( Gf[0][0] - 1. ), std::fabs( Gf[1][0] ), std::fabs( Ga[0][0] - 1. ), std::fabs( Ga[1][0] ) } ): 1.;
      check( std::string( "(G) pointwise " ) + ( which == 0? "DAESLV": "ODESLV" ) + ": x(0.5+) = x(0.5-) + g, J = 1.3 and dJ/dg = (1, 0), forward and adjoint"
             + ( ok? std::string(): " -- " + w ), ok && err < 1e-8 );
    }
  }
  std::cout << "==========  DAESLV_moltime: " << ( nfail? std::to_string( nfail ) + " FAILED": std::string( "ALL PASS" ) ) << "  ==========\n";
  return nfail? 1: 0;
}
