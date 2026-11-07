// OCFE_byname.cpp -- OCFESLV's by-name solve, the form ODESLV_CVODES and DAESLV_IDAS share:
//   oc.solve({ {k0, {0.5}}, {c_ic, []( FFModel::t_Coord const& c ){ ... }} });  oc.val_function();
// Gates: (A) it reproduces the in-place path (init + set_input_values + sample_input + solve) exactly, under the
// monolithic and the marching solve, outputs included (a terminal spatial integral, an outlet time integral);
// (B) a coordinate generator and the equivalent flat list in input_dofs() order give the same result;
// (C) an incomplete or inconsistent list is refused, naming the entry.
#include <iostream>
#include <iomanip>
#include <cmath>
#include "ocfeslv.hpp"
using namespace mc;
static int npass = 0, nfail = 0;
static void check( bool ok, std::string const& what ){ std::cout << "  " << std::left << std::setw(68) << what << ( ok? "PASS": "FAIL" ) << "\n"; ( ok? npass: nfail )++; }

int main()
{
  double const D = 0.05, K0 = 0.5, T_end = 1.;
  auto f_ic = []( double z ){ return 0.5 + 0.5*std::cos( M_PI*z ); };       // c_ic(0) = 1 = the inflow value
  for( int march = 0; march < 2; ++march ){
    std::cout << ( march? "MARCHING": "MONOLITHIC" ) << "\n";
    FFGraph G;  OCFESLV oc( &G );  FFPartial OpP;  FFIntegral OpI;
    FFVar t = G.add_var( "t" ), z = G.add_var( "z" ), c = G.add_var( "c" ), k0 = G.add_var( "k0" ), c_ic = G.add_var( "c_ic" );
    oc.add_domain( t, FFDom( 0., T_end, 3, FFDom::LGR, 3 ) );
    oc.add_domain( z, FFDom( 0., 1., 3, FFDom::CGL, 5 ) );
    oc.add_state( c, {t,z} );
    oc.add_input( k0, {} );                                                 // no reference: must be given
    oc.add_input( c_ic, {z} );
    oc.update_ref( c_ic, [&]( OCFESLV::t_Coord const& ){ return 0.; } );
    oc.set_evolution_domain( t );
    int const T_NO_LB = FFDom::ALL - FFDom::LB, Z_INT = FFDom::ALL - FFDom::LB - FFDom::UB;
    oc.add_equation( OpP( c, t ) - D*OpP( c, {z,2} ) + k0*c, {t,z}, {T_NO_LB, Z_INT}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INTERIOR, 0 ) );
    oc.add_equation( c - c_ic,       {t,z}, {FFDom::LB, FFDom::ALL}, OCFESLV::EqnOptions( OCFESLV::EqnRole::INITIAL,  0 ) );
    oc.add_equation( c - 1.,         {t,z}, {T_NO_LB, FFDom::LB},    OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
    oc.add_equation( OpP( c, z ),    {t,z}, {T_NO_LB, FFDom::UB},    OCFESLV::EqnOptions( OCFESLV::EqnRole::BOUNDARY, 0 ) );
    oc.add_output( OpI( c, z ), {t}, {T_end} );                             // fct 0: terminal inventory
    oc.add_output( OpI( c, t ), {z}, {1.0} );                               // fct 1: outlet time integral
    oc.options.REDUCE.ORDER = OCFESLV::Options::RED_FULL;
    oc.options.DISPLAY_LEVEL = 0;
    if( march ) oc.options.SOLVE.WARMSTART = OCFESLV::Options::BROADCAST_IC;
    else        oc.options.SOLVE.MARCHING  = false;
    if( !oc.setup() ){ check( false, "setup" ); continue; }

    // (A) the in-place path
    std::vector<double> xv, inp;
    oc.init( xv, inp, nullptr );
    oc.set_input_values( k0, { K0 }, inp.data() );
    oc.sample_input( c_ic, [&]( OCFESLV::t_Coord const& cr ){ return f_ic( cr.at( z ) ); }, inp.data() );
    auto const rep = oc.solve( xv.data(), inp.data(), nullptr );
    std::vector<double> const ref = oc.val_functions();
    check( rep.converged && ref.size() == 2, "in-place solve converges, 2 outputs" );

    // (A) by name, the coordinate generator
    auto const st = oc.solve({ { k0, { K0 } }, { c_ic, FFModel::t_Fun( [&]( FFModel::t_Coord const& cr ){ return f_ic( cr.at( z ) ); } ) } });
    double dA = 0.; for( size_t i = 0; i < ref.size() && i < oc.val_function().size(); ++i ) dA = std::max( dA, std::fabs( oc.val_function()[i] - ref[i] ) );
    double dX = 0.; for( size_t i = 0; i < xv.size() && i < oc.val_state().size(); ++i ) dX = std::max( dX, std::fabs( oc.val_state()[i] - xv[i] ) );
    std::cout << std::scientific << std::setprecision(3) << "    Inv(T) = " << ref[0] << "   Eff(z=1) = " << ref[1]
              << "   |by name - in place|: outputs " << dA << ", state " << dX << "\n";
    check( st == OCFESLV::NORMAL && dA == 0. && dX == 0., "(A) by name = in place, bit for bit (outputs and state)" );

    // (B) the flat list in input_dofs() order, from input_coordinates
    std::vector<double> flat;  for( auto const& cr : oc.input_coordinates( c_ic ) ) flat.push_back( f_ic( cr.at( z ) ) );
    oc.solve({ { k0, { K0 } }, { c_ic, flat } });
    double dB = 0.; for( size_t i = 0; i < ref.size(); ++i ) dB = std::max( dB, std::fabs( oc.val_function()[i] - ref[i] ) );
    check( dB == 0. && flat.size() == oc.input_dofs( c_ic ).size(), "(B) the flat list = the coordinate generator" );

    // (C) refusals name the entry
    auto refused = [&]( std::vector<FFModel::InputVal> const& v, std::string const& must ){
      try{ oc.solve( v ); }catch( std::invalid_argument const& e ){ return std::string( e.what() ).find( must ) != std::string::npos; }
      return false; };
    check( refused( { { k0, { K0 } } }, "missing: c_ic" ),                        "(C) a missing input is refused, named" );
    check( refused( { { k0, { K0 } }, { c_ic, { 1., 2. } } }, "c_ic has " ),        "(C) a wrong count is refused, named" );
    check( refused( { { k0, { K0 } }, { k0, { K0 } }, { c_ic, flat } }, "k0 is given more than once" ), "(C) a duplicate is refused, named" );
    // (D) SENSITIVITIES by name (phase 7): forward and adjoint, one row per declared DOF (k0's one, then c_ic's in
    //     input_dofs order), one column per output -- = the registry path's sens_jacobian, = central differences
    {
      std::vector<FFModel::InputVal> const in{ { k0, { K0 } }, { c_ic, flat } };
      auto const sf = oc.solve_fsens( in );
      auto const Gf = oc.val_function_gradient();
      auto const Ff = oc.val_function();
      size_t const nd = 1 + oc.input_dofs( c_ic ).size();
      check( sf == OCFESLV::NORMAL && Gf.size() == nd && Gf.size() && Gf[0].size() == 2,
             "(D) solve_fsens by name: one row per declared DOF, one column per output" );
      auto const sa = oc.solve_asens( in );
      auto const Ga = oc.val_function_gradient();
      double dfa = 0.;
      for( size_t r = 0; r < Gf.size() && r < Ga.size(); ++r ) for( size_t f = 0; f < 2; ++f ) dfa = std::max( dfa, std::fabs( Gf[r][f] - Ga[r][f] ) );
      check( sa == OCFESLV::NORMAL && Ga.size() == Gf.size() && dfa < 1e-9, "(D) solve_asens by name = solve_fsens, to 1e-9" );
      // the registry path: the same Jacobian, laid out over the control vector
      oc.register_control( k0 );  oc.register_control( c_ic );
      std::vector<double> xs, is;  oc.init( xs, is, nullptr );
      oc.set_input_values( k0, { K0 }, is.data() );  oc.set_input_values( c_ic, flat, is.data() );
      bool const okr = oc.solve_fsens( xs.data(), is.data(), nullptr );
      auto const& J = oc.sens_jacobian();  size_t const ncd = oc.n_control_dof();
      double dreg = 0.;
      if( okr && ncd == nd ){
        auto const bk = oc.control_block( k0 ), bc = oc.control_block( c_ic );
        for( size_t f = 0; f < 2; ++f ){
          dreg = std::max( dreg, std::fabs( Gf[0][f] - J[ f*ncd + bk.offset ] ) );
          for( size_t k = 0; k < bc.ndof; ++k ) dreg = std::max( dreg, std::fabs( Gf[1+k][f] - J[ f*ncd + bc.offset + k ] ) );
        }
      }
      check( okr && ncd == nd && dreg == 0., "(D) by name = the registry path's sens_jacobian, bit for bit" );
      // central differences on k0 and on one c_ic DOF, against the forward rows
      double const h = 1e-5;  double dfd = 0.;
      auto J_at = [&]( double k, size_t kdof, double dk ){
        std::vector<double> fl = flat;  if( kdof < fl.size() ) fl[kdof] += dk;
        oc.solve( { { k0, { k } }, { c_ic, fl } } );  return oc.val_function(); };
      { auto fp = J_at( K0 + h, flat.size(), 0. ), fm = J_at( K0 - h, flat.size(), 0. );
        for( size_t f = 0; f < 2; ++f ) dfd = std::max( dfd, std::fabs( ( fp[f] - fm[f] ) / ( 2*h ) - Gf[0][f] ) ); }
      { size_t const kd = flat.size() / 2;
        auto fp = J_at( K0, kd, h ), fm = J_at( K0, kd, -h );
        for( size_t f = 0; f < 2; ++f ) dfd = std::max( dfd, std::fabs( ( fp[f] - fm[f] ) / ( 2*h ) - Gf[1+kd][f] ) ); }
      std::cout << "    d(outputs)/dk0: " << Gf[0][0] << ", " << Gf[0][1] << "   max |FD - fwd| " << dfd << "\n";
      check( dfd < 1e-6, "(D) forward rows = central differences (k0, a c_ic DOF), to 1e-6" );
      oc.clear_controls();
    }
  }
  std::cout << "\n==========  OCFE_byname: " << npass << " passed, " << nfail << " failed" << ( nfail? "": " -- ALL PASS" ) << "  ==========\n";
  return nfail? 1: 0;
}
