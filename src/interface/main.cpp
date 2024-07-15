#include <pybind11/pybind11.h>

namespace py = pybind11;

void mc_odeslvs( py::module_ & );
void mc_ffunc( py::module_ & );

PYBIND11_MODULE( cronos, m )
{

  m.doc() = "Python interface of library CRONOS";

  mc_odeslvs( m );
  mc_ffunc( m );

}

