// Checks of mc::FFBLOCK (ffblock.hpp): block lower-triangular and bordered lower-triangular forms
#include <random>
#include <numeric>
#include <functional>
#include "ffblock.hpp"

static int nfail = 0, nok = 0;
#define CHECK(c,m) do{ if(c) ++nok; else { ++nfail; std::cout << "FAIL " << m << std::endl; } }while(0)

static bool has( std::vector<int> const& r, int j ){ return std::find( r.begin(), r.end(), j ) != r.end(); }

// independent maximum matching (augmenting paths), rows -> columns
static int maxmatch( unsigned m, unsigned n, std::vector<std::vector<int>> const& rows )
{
  std::vector<int> mcol( n, -1 );
  std::function<bool(int,std::vector<char>&)> aug = [&]( int i, std::vector<char>& seen ){
    for( int j : rows[i] ){ if( seen[j] ) continue; seen[j] = 1;
      if( mcol[j] < 0 || aug( mcol[j], seen ) ){ mcol[j] = i; return true; } }
    return false; };
  int c = 0; for( unsigned i=0; i<m; i++ ){ std::vector<char> s( n, 0 ); if( aug( i, s ) ) c++; }
  return c;
}

// independent strongly connected components (Kosaraju)
static int nscc( unsigned n, std::vector<std::vector<int>> const& adj, std::vector<int>& comp )
{
  std::vector<std::vector<int>> radj( n );
  for( unsigned k=0; k<n; k++ ) for( int l : adj[k] ) radj[l].push_back( k );
  std::vector<char> seen( n, 0 ); std::vector<int> order;
  std::function<void(int)> d1 = [&]( int k ){ seen[k] = 1; for( int l : adj[k] ) if( !seen[l] ) d1( l ); order.push_back( k ); };
  for( unsigned k=0; k<n; k++ ) if( !seen[k] ) d1( k );
  comp.assign( n, -1 ); int c = 0;
  std::function<void(int)> d2 = [&]( int k ){ comp[k] = c; for( int l : radj[k] ) if( comp[l] < 0 ) d2( l ); };
  for( auto it=order.rbegin(); it!=order.rend(); ++it ) if( comp[*it] < 0 ){ d2( *it ); c++; }
  return c;
}

static bool isperm( int const* p, unsigned n ){
  std::vector<int> a( p, p+n ); std::sort( a.begin(), a.end() );
  for( unsigned i=0; i<n; i++ ) if( a[i] != (int)i+1 ) return false; return true; }

static void check_btf( unsigned n, std::vector<std::vector<int>> const& rows, std::string const& tag )
{
  std::vector<int> IPERM( n ), IOR( n ), IB( n ); int NB = -1;
  bool const ok = mc::FFBLOCK::block_triangular( n, rows, IPERM.data(), IOR.data(), IB.data(), NB );
  CHECK( ok == ( maxmatch( n, n, rows ) == (int)n ), tag+" singularity" );
  if( !ok ) return;
  CHECK( isperm( IOR.data(), n ) && isperm( IPERM.data(), n ), tag+" permutations" );
  std::vector<int> rowat( n ), colpos( n ), blk( n );
  for( unsigned k=0; k<n; k++ ){ rowat[k] = IPERM[IOR[k]-1]-1; colpos[IOR[k]-1] = k; }
  for( int b=0; b<NB; b++ ){ int e = b+1<NB? IB[b+1]-1: n; for( int k=IB[b]-1; k<e; k++ ) blk[k] = b; }
  bool diag = true, lower = true; std::vector<std::vector<int>> adj( n );
  for( unsigned k=0; k<n; k++ ){
    diag = diag && has( rows[rowat[k]], IOR[k]-1 );
    for( int j : rows[rowat[k]] ){ int l = colpos[j]; if( blk[l] > blk[k] ) lower = false; if( l != (int)k ) adj[k].push_back( l ); }
  }
  std::vector<int> comp; int nc = nscc( n, adj, comp ); bool same = nc == NB;
  for( unsigned k=0; same && k<n; k++ ) for( unsigned l=0; same && l<n; l++ ) same = ( comp[k]==comp[l] ) == ( blk[k]==blk[l] );
  CHECK( diag && lower && same, tag+" zero-free diagonal, block lower-triangular, blocks = SCCs" );
}

#ifndef CRONOS__WITH_HSL
// exact minimum number of pairs to tear so that the matched dependency graph is acyclic
static int min_tear( std::vector<std::vector<int>> const& dep )
{
  unsigned const np = dep.size();
  for( unsigned t=0; t<=np; t++ ){
    std::vector<int> sel( np, 0 ); std::fill( sel.end()-t, sel.end(), 1 );
    do{
      std::vector<int> indeg( np, 0 ); unsigned cnt = 0;
      for( unsigned k=0; k<np; k++ ) if( !sel[k] ) for( int l : dep[k] ) if( !sel[l] ) indeg[k]++;
      std::vector<int> st; for( unsigned k=0; k<np; k++ ) if( !sel[k] && !indeg[k] ) st.push_back( k );
      std::vector<std::vector<int>> user( np ); for( unsigned k=0; k<np; k++ ) for( int l : dep[k] ) user[l].push_back( k );
      while( !st.empty() ){ int k = st.back(); st.pop_back(); cnt++; for( int u : user[k] ) if( !sel[u] && !--indeg[u] ) st.push_back( u ); }
      if( cnt == np - t ) return t;
    } while( std::next_permutation( sel.begin(), sel.end() ) );
  }
  return np;
}
#endif

static void check_bordered( unsigned m, unsigned n, std::vector<std::vector<int>> const& rows, std::string const& tag,
                            long* sumb=nullptr, long* sumopt=nullptr )
{
  std::vector<int> IP( m ), IQ( n ); int NROW, NCOL, NBROW, NBCOL;
  bool const ok = mc::FFBLOCK::bordered_triangular( m, n, rows, IP.data(), IQ.data(), NROW, NCOL, NBROW, NBCOL );
  unsigned nnz = 0; for( auto& r : rows ) nnz += r.size();
  CHECK( ok == ( nnz > 0 ), tag+" bordered: status" );
  if( !ok ) return;
  CHECK( isperm( IP.data(), m ) && isperm( IQ.data(), n ), tag+" bordered: permutations" );
  std::vector<int> colpos( n ); for( unsigned j=0; j<n; j++ ) colpos[IQ[j]-1] = j;
  std::vector<char> colnz( n, 0 ); for( auto& r : rows ) for( int j : r ) colnz[j] = 1;
  bool nulls = true;
  for( unsigned i=0; i<m; i++ ) nulls = nulls && ( (int)i < NROW ) == !rows[IP[i]-1].empty();
  for( unsigned j=0; j<n; j++ ) nulls = nulls && ( (int)j < NCOL ) == (bool)colnz[IQ[j]-1];
  CHECK( nulls, tag+" bordered: null rows and columns last" );
  int const LR = NROW-NBROW, LC = NCOL-NBCOL;
#ifndef CRONOS__WITH_HSL
  // BTF: maximum transversal and tearing -- 1x1 diagonal blocks
  bool lower = LR >= LC, diag = true;
  for( int i=0; i<LR; i++ ){
    auto const& r = rows[IP[i]-1];
    for( int j : r ){ int c = colpos[j]; if( c < LC && c > i ) lower = false; }
    if( i < LC ) diag = diag && has( r, IQ[i]-1 );
  }
  CHECK( lower && diag, tag+" bordered: leading block lower-triangular with zero-free diagonal" );
#else
  // HSL MC33 (P5): block lower-triangular with full diagonal blocks -- the blocks are the runs tied
  // by the entries above the diagonal; each must have a zero-free diagonal
  std::vector<int> reach( std::max( LR, LC ) );
  std::iota( reach.begin(), reach.end(), 0 );
  for( int i=0; i<LR; i++ ) for( int j : rows[IP[i]-1] ){ int c = colpos[j]; if( c < LC && c > i ) reach[i] = std::max( reach[i], c ); }
  bool diag = true;
  for( int s=0, e=0; s<std::min( LR, LC ); s=e+1 ){
    e = reach[s]; for( int k=s; k<=e && k<(int)reach.size(); k++ ) e = std::max( e, reach[k] );
    for( int k=s; k<=e && k<std::min( LR, LC ); k++ ) diag = diag && has( rows[IP[k]-1], IQ[k]-1 );
  }
  CHECK( diag, tag+" bordered: leading block block-lower-triangular with zero-free diagonal" );
#endif
  // the leading square block is as large as a maximum matching allows, less the torn pairs
  CHECK( LC + NBCOL >= 0 && LC <= maxmatch( m, n, rows ), tag+" bordered: leading block within the matching" );
#ifndef CRONOS__WITH_HSL
  if( sumb && m == n && maxmatch( m, n, rows ) == (int)n && n <= 9 ){
    // exact minimum border for the same matching (pairs = leading + torn)
    std::vector<int> pairrow, paircol, colpair( n, -1 );
    for( int k=0; k<LC; k++ ){ pairrow.push_back( IP[k]-1 ); paircol.push_back( IQ[k]-1 ); }
    for( int k=LR; k<NROW; k++ ){ pairrow.push_back( IP[k]-1 ); paircol.push_back( IQ[LC+k-LR]-1 ); }
    for( unsigned k=0; k<pairrow.size(); k++ ) colpair[paircol[k]] = k;
    std::vector<std::vector<int>> dep( pairrow.size() );
    for( unsigned k=0; k<pairrow.size(); k++ ) for( int j : rows[pairrow[k]] ) if( colpair[j] >= 0 && colpair[j] != (int)k ) dep[k].push_back( colpair[j] );
    *sumb += NBCOL; *sumopt += min_tear( dep );
  }
#endif
}

int main()
{
  std::mt19937 gen( 7 );
  long sumb = 0, sumopt = 0; int nq = 0;
  for( int t=0; t<3000; t++ ){
    bool const square = t%3;
    unsigned const m = 1 + gen()%40, n = square? m: 1 + gen()%40;
    double const dens = std::uniform_real_distribution<>( 0., 0.35 )( gen );
    std::vector<int> perm( n ); std::iota( perm.begin(), perm.end(), 0 ); std::shuffle( perm.begin(), perm.end(), gen );
    std::vector<std::vector<int>> rows( m );
    for( unsigned i=0; i<m; i++ ){
      if( i < n && gen()%8 ) rows[i].push_back( perm[i] );
      for( unsigned j=0; j<n; j++ ) if( std::uniform_real_distribution<>()( gen ) < dens ) rows[i].push_back( j );
      if( gen()%25 == 0 ) rows[i].clear(); // null row
      std::sort( rows[i].begin(), rows[i].end() ); rows[i].erase( std::unique( rows[i].begin(), rows[i].end() ), rows[i].end() );
    }
    std::string tag = "#"+std::to_string(t)+" "+std::to_string(m)+"x"+std::to_string(n);
    if( square ) check_btf( n, rows, tag );
    if( square && n <= 9 ) nq++;
    check_bordered( m, n, rows, tag, &sumb, &sumopt );
  }
#ifndef CRONOS__WITH_HSL
  std::cout << "border size of the greedy tearing vs exact minimum (same matching, " << nq << " square systems n<=9): "
            << sumb << " vs " << sumopt << std::endl;
#endif

  // DAG interface
  mc::FFGraph DAG; mc::FFVar P[2], X[4];
  for( auto& p : P ) p.set( &DAG ); for( auto& x : X ) x.set( &DAG );
  std::vector<mc::FFVar> sys = { X[0] - X[3]*P[0], X[1] - X[0]*X[0], X[2] - X[1] - X[0], X[3] - X[2]*P[1] }; // one cycle
  mc::FFBLOCK blk( &DAG );
  int IPERM[4], IOR[4], IB[4], NB, IP[4], IQ[4], NROW, NCOL, NBROW, NBCOL;
  CHECK( blk.block_triangular( 4, sys.data(), X, IPERM, IOR, IB, NB, true ) && NB == 1, "cyclic system: 1 block" );
  CHECK( blk.bordered_triangular( 4, sys.data(), 4, X, IP, IQ, NROW, NCOL, NBROW, NBCOL, true ) && NBROW == 1 && NBCOL == 1, "cyclic system: one tear" );
  std::vector<mc::FFVar> sys2 = { X[0] + X[1] - P[0], X[1]*X[2] - P[1], X[0] - X[2] + X[3] }; // 3 equations, 4 variables
  CHECK( blk.bordered_triangular( 3, sys2.data(), 4, X, IP, IQ, NROW, NCOL, NBROW, NBCOL, true ), "underdetermined system" );

  std::cout << mc::FFBLOCK::library() << ": " << nok << " ok, " << nfail << " failed" << std::endl;
  return nfail? 1: 0;
}
