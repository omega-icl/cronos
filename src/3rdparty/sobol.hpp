/*
                                                  SOBOL
                                     The Sobol Quasirandom Sequence
                      (http://people.sc.fsu.edu/~jburkardt/cpp_src/sobol/sobol.html)

SOBOL is a C++ library which computes elements of the Sobol quasirandom sequence, by Bennett Fox.

A quasirandom or low discrepancy sequence, such as the Faure, Halton, Hammersley, Niederreiter or Sobol sequences, is "less random" than a pseudorandom number sequence, but more useful for such tasks as approximation of integrals in higher dimensions, and in global optimization. This is because low discrepancy sequences tend to sample space "more uniformly" than random numbers. Algorithms that use such sequences may have superior convergence.

SOBOL is an adapation of the INSOBL and GOSOBL routines in ACM TOMS Algorithm 647 and ACM TOMS Algorithm 659. The original code can only compute the "next" element of the sequence. The revised code allows the user to specify the index of any desired element.

A remark by Joe and Kuo shows how to extend the algorithm from the original maximum spatial dimension of 40 up to a maximum spatial dimension of 1111. The FORTRAN90 and C++ versions of this program have been updated in this way. In particular, the extra data in the C++ version of the program was kindly formatted and supplied by Steffan Berridge.

The routines with a prefix of I8_ use 64 bit integers, and use the long int to get this. On some systems, a long int is simply 32 bits. In that case, try using the long long int datatype instead.

The original, true, correct versions of ACM TOMS Algorithm 647 and ACM TOMS Algorithm 659 are available in the TOMS subdirectory of the NETLIB web site. The version displayed here has been converted to FORTRAN90, and other internal changes have been made to suit me.
Licensing:

The computer code and data files described and made available on this web page are distributed under the GNU LGPL license.
Languages:

SOBOL is available in a C++ version and a FORTRAN90 version and a MATLAB version and a Python version 
*/
int i4_bit_hi1 ( int n );
int i4_bit_lo0 ( int n );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void i4_sobol ( int dim_num, int *seed, float quasi[ ] );
float *i4_sobol_generate ( int m, int n, int skip );
int i4_uniform ( int b, int c, int *seed );

int i8_bit_hi1 ( long long int n );
int i8_bit_lo0 ( long long int n );
long long int i8_max ( long long int i1, long long int i2 );
long long int i8_min ( long long int i1, long long int i2 );
void i8_sobol ( int dim_num, long long int *seed, double quasi[ ] );
double *i8_sobol_generate ( int m, int n, int skip );
long long int i8_uniform ( long long int b, long long int c, int *seed );

float r4_abs ( float x );
int r4_nint ( float x );
float r4_uniform_01 ( int *seed );

double r8_abs ( double x );
int r8_nint ( double x );
double r8_uniform_01 ( int *seed );

void r8mat_write ( std::string output_filename, int m, int n, double table[] );

int tau_sobol ( int dim_num );
void timestamp ( );

