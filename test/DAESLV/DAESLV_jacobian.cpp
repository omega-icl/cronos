// DAESLV_jacobian -- the analytic Jacobian dF/dy + c_j dF/dy', checked against the residual itself.
// =====================================================================================================
//  F  Backend agreement is weak evidence -- a wrong Jacobian usually still converges, just more slowly.
//     This differences the residual on a NONLINEAR model, where every entry is also hand-checkable.
#include <cmath>
#include <iomanip>
#include <iostream>
#include "daeslv_idas.hpp"
#include "odeslv_cvodes.hpp"

using namespace mc;

//! Exposes the two callbacks so the ANALYTIC Jacobian can be checked against finite differences of the
//! residual itself -- agreement between solver backends is weak evidence, since a wrong Jacobian usually
//! still converges, just more slowly.
struct JacProbe : public DAESLV_IDAS {
  int res( double t, N_Vector y, N_Vector yp, N_Vector r ){ return IDRES__( t, y, yp, r, nullptr ); }
  int jac( double t, double cj, N_Vector y, N_Vector yp, SUNMatrix J )
    { return IDJAC__( t, cj, y, yp, nullptr, J, nullptr, nullptr, nullptr, nullptr ); }
};

static int n_pass = 0, n_fail = 0;
static void check( std::string const& w, double got, double want, double tol )
{
  double const d = std::fabs( got - want );
  ( d <= tol )? ++n_pass: ++n_fail;
  std::cout << "    " << std::left << std::setw( 40 ) << w << std::right
            << " got=" << std::scientific << std::setprecision(10) << got
            << " want=" << want << " |d|=" << std::setprecision(2) << d
            << ( d <= tol? "  PASS\n": "  FAIL\n" );
}
static double const T_END = 1.0;

int main()
{
  // ---------------------------------------------------------------- F
  { std::cout << "\n--- F. analytic Jacobian vs finite differences of the residual (NONLINEAR model)\n";
    FFGraph G; JacProbe M; M.set_dag( &G );
    FFVar t=G.add_var("t"), y=G.add_var("y(t)"), z=G.add_var("z(t)");
    FFPartial OpP; int const TI = FFDom::ALL - FFDom::LB;
    M.add_domain( t, FFDom( 0., 1., 2, FFDom::LGR, 4 ) );
    M.add_state( y, {t} ); M.add_state( z, {t} );
    M.set_evolution_domain( t ); M.update_ref( y, 1. ); M.update_ref( z, 1. );
    auto io=FFModel::EqnOptions(FFModel::EqnRole::INTERIOR,0), ii=FFModel::EqnOptions(FFModel::EqnRole::INITIAL,0);
    M.add_equation( OpP(y,t) + y*z, {t}, {TI},         io );   // nonlinear in BOTH states
    M.add_equation( z - y*y,        {t}, {FFDom::ALL}, io );
    M.add_equation( y - 1.,         {t}, {FFDom::LB},  ii );
    M.options.DISPLAY_LEVEL = 0;
    if( !M.setup() ){ std::cout << "    setup failed: " << M.extract_error() << "\n"; ++n_fail; }
    else{
      SUNContext ctx; SUNContext_Create( SUN_COMM_NULL, &ctx );
      size_t const n = M.id().size();
      N_Vector yy=N_VNew_Serial(n,ctx), yp=N_VNew_Serial(n,ctx), r0=N_VNew_Serial(n,ctx), r1=N_VNew_Serial(n,ctx);
      SUNMatrix J = SUNDenseMatrix(n,n,ctx);
      double const tt = 0.37, cj = 2.75;                      // arbitrary, away from any special point
      N_VGetArrayPointer(yy)[0]=0.83; N_VGetArrayPointer(yy)[1]=0.61;
      N_VGetArrayPointer(yp)[0]=-0.42; N_VGetArrayPointer(yp)[1]=0.19;
      M.jac( tt, cj, yy, yp, J );
      M.res( tt, yy, yp, r0 );
      double worst = 0.;
      double const h = 1e-7;
      for( size_t j = 0; j < n; ++j ){
        // dF/dy_j
        double const sy = N_VGetArrayPointer(yy)[j];
        N_VGetArrayPointer(yy)[j] = sy + h;  M.res( tt, yy, yp, r1 );  N_VGetArrayPointer(yy)[j] = sy;
        std::vector<double> fd(n);
        for( size_t k = 0; k < n; ++k ) fd[k] = ( N_VGetArrayPointer(r1)[k] - N_VGetArrayPointer(r0)[k] ) / h;
        // + c_j dF/dyp_j
        double const sp = N_VGetArrayPointer(yp)[j];
        N_VGetArrayPointer(yp)[j] = sp + h;  M.res( tt, yy, yp, r1 );  N_VGetArrayPointer(yp)[j] = sp;
        for( size_t k = 0; k < n; ++k )
          fd[k] += cj * ( N_VGetArrayPointer(r1)[k] - N_VGetArrayPointer(r0)[k] ) / h;
        for( size_t k = 0; k < n; ++k ){
          double const d = std::fabs( SM_ELEMENT_D(J,k,j) - fd[k] );
          if( d > worst ) worst = d;
          std::cout << "      J(" << k << "," << j << ") analytic=" << std::setw(12) << std::setprecision(6)
                    << std::fixed << SM_ELEMENT_D(J,k,j) << "   fd=" << std::setw(12) << fd[k] << "\n";
        }
      }
      std::cout << std::scientific;
      check( "analytic Jacobian matches finite differences", worst, 0., 1e-5 );
      SUNMatDestroy(J); N_VDestroy(yy); N_VDestroy(yp); N_VDestroy(r0); N_VDestroy(r1); SUNContext_Free(&ctx);
    } }


  std::cout << "\n==========  DAESLV_jacobian: " << n_pass << " passed, " << n_fail << " failed -- "
            << ( n_fail? "FAILURES": "ALL PASS" ) << "  ==========\n\n";
  return n_fail? 1: 0;
}
