// ODESLV_revision.cpp -- the PROVENANCE of an ODESLV sweep: prints the revisions of the headers this binary was built
// with, once per sweep log -- the model layer (ffmodel.hpp) and the ODESLV headers (odeslv_base.hpp and the three that
// change with it).  Before 2026-10-07 the ODESLV headers carried no revision, so a log could not show which ran.
#include <cstdio>
#include <cstring>
#include <string>
#include "ffode.hpp"

static std::string rev( char const* s ){                 // "odeslv  rev1  2026-10-07" -> "rev1"
  std::string t( s ); auto p = t.find( "rev" );
  if( p == std::string::npos ) return std::string();
  auto q = t.find( ' ', p );
  return t.substr( p, q == std::string::npos? q: q - p );
}

int main(){
  std::printf( "CRONOS header revisions\n" );
  std::printf( "  model layer (ffmodel.hpp)     : %s\n", mc::FFModel::revision() );
  std::printf( "  ODE solver  (odeslv*.hpp)     : %s\n", mc::ODESLVS_CVODES::revision() );
  bool const ok = !rev( mc::FFModel::revision() ).empty() && !rev( mc::ODESLVS_CVODES::revision() ).empty()
               && std::strcmp( mc::FFModel::revision(), mc::ODESLVS_CVODES::revision() ) != 0;
  std::printf( "  %s\n", ok? "PASS  both revisions present": "FAIL  a revision is missing" );
  return ok? 0: 1;
}
