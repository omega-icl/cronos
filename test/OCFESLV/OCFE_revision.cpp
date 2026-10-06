// OCFE_revision.cpp -- the PROVENANCE of a sweep: prints the revisions of the CRONOS headers this binary was built
// with, once per sweep log.  Since 2026-10-04 the OCFESLV setup banner prints them only at DISPLAY_LEVEL >= 2;
// revision() is the query (FFModel::revision(), OCFESLV::revision()).  Both are printed because the model layer
// and the solver are separate headers with their own revision numbers, and a build can combine any two.
#include <cstdio>
#include <cstring>
#include <string>
#include "ocfeslv.hpp"

static std::string rev( char const* s ){                 // "ocfeslv  rev354  2026-10-04" -> "rev354"
  std::string t( s ); auto p = t.find( "rev" );
  if( p == std::string::npos ) return std::string();
  auto q = t.find( ' ', p );
  return t.substr( p, q == std::string::npos? q: q - p );
}

int main(){
  std::printf( "CRONOS header revisions\n" );
  std::printf( "  model layer (ffmodel.hpp) : %s\n", mc::FFModel::revision() );
  std::printf( "  solver      (ocfeslv.hpp) : %s\n", mc::OCFESLV::revision() );
  bool const ok = std::strlen( mc::FFModel::revision() ) && std::strlen( mc::OCFESLV::revision() )
               && !rev( mc::FFModel::revision() ).empty() && !rev( mc::OCFESLV::revision() ).empty();
  std::printf( "  %s\n", ok? "PASS  both revisions present": "FAIL  a revision is missing" );
  return ok? 0: 1;
}
