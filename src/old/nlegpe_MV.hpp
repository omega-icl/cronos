#ifndef MC__NLEGPE_H
#define MC__NLEGPE_H

#include <iostream>
#include <list>
#include "base_nle.hpp"
#include "setinv.hpp"
#include "cmodel.hpp"
#include "ellipsoid.hpp"
#include "gurobi_c++.h"

#undef  MC__NLEGPE_DEBUG


namespace mc
{

template <typename T>
class NLEGPE: protected SetInv<T>, public BASE_NLE
{
 public:

  typedef mc::CModel<T> CM;
  typedef mc::CVar<T> CV;

  //Constructor
  NLEGPE(): _nym(0), _pGRB(0), _pCM(0)
    {};

  // Default Destructor
  ~NLEGPE()
    { delete _pGRB; delete _pCM; }

  //! @brief NLEGPE options
  struct Options
  {
    //! @brief Constructor
    Options():
      SETINV(typename SetInv<T>::Options()), CMODEL(typename CM::Options()),
      OUTPUT_BOUND(CM), CM_ORDER(2), OUTRED_MAX(10), OUTRED_THRES(2e-2),
      OUTRED_TOL(1e-9), INRED_MAX(10), INRED_THRES(1e-2), ROTATION_USE(false)
      {}
    //! @brief Display
    void display
      ( std::ostream&out ) const;
    //! @brief output bounder
    enum BOUNDER{
      IA=0,	//!< Interval analysis
      CM	//!< Chebyshev model
    };
    //! @brief Options of SetInv class
    typename SetInv<T>::Options SETINV;
    //! @brief Options of CModel class
    typename CM::Options CMODEL;
    //! @brief output bounding strategy
    BOUNDER OUTPUT_BOUND;
    //! @brief Order of Chebyshev model for output bounding
    unsigned CM_ORDER;
    //! @brief Maximum number of domain outer-reduction
    unsigned OUTRED_MAX;
    //! @brief Threshold for variable domain outer-reduction (minimum domain reduction ratio)
    double OUTRED_THRES;
    //! @brief LP solver tolerance for domain outer-reduction
    double OUTRED_TOL;
   //! @brief Maximum number of domain inner-reduction
    unsigned INRED_MAX;
    //! @brief Threshold for variable domain inner-reduction (minimum domain reduction ratio)
    double INRED_THRES;
     //! @brief Whether to rotate variable domains
    bool ROTATION_USE;
  } options;

  //! @brief NLEGPE exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for NLEGPE exception handling
    enum TYPE{
      INTERN=-33	//!< Internal error
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Inline function returning the error flag
    int ierr(){ return _ierr; }
    std::string what(){
      switch( _ierr ){
      case INTERN: default:
        return "Internal error";
      }
    }
  private:
    TYPE _ierr;
  };

  //! @brief Structure holding current statistics
  struct Stats{
    double BOUND;
    double OUTRED;
    double INRED;
    double ROTATE;
    double SETINV;
  } stats;

  //! @brief Structure to hold measurement data
  struct Data{
    Data( double _yl, double _yu, unsigned _iy=0, double _t=0. ):
      t(_t), iy(_iy), yl(_yl), yu(_yu)
      {}
    double   t;
    unsigned iy;
    double   yl;
    double   yu;
  };

  //! @brief Solve NLEGPE problem -- return value is pair of inner-approximation and boundary volumes
  double solve
    ( const T*P, const std::list<Data>&Ym, std::ostream&os=std::cout );

  //! @brief Append all open nodes and inner nodes to, respectively, <a>os_open</a> and <a>os_inner</a> - Number of significant digits is set via <a>DPREC</a> (default=6)
  void output_nodes
    ( std::ostream&os_open, const bool PROJ=false, const unsigned int DPREC=6 ) const;

  //! @brief Return the current stack of (open) nodes
  const typename SetInv<T>::t_Nodes& open_nodes
    () const;

  //! @brief Cluster boxes in set-boundary approximation - on return, list of box enclosures for clusters
  std::list<T*> clusters
    ( const double rtol=1e-4, const double atol=1e-8, std::ostream&os=std::cout ) const;

  //! @brief Append all clusters to <a>os_clus</a> - Number of significant digits is set via <a>DPREC</a> (default=6)
  void output_clusters
    ( std::ostream&os_clus, const double RTOL=1e-4, const double ATOL=1e-8,
      const unsigned DPREC=6 ) const;
  
  //! @brief Append all clusters to <a>os_clus</a> - Number of significant digits is set via <a>DPREC</a> (default=6)    
  void output_node
	( const SetInvNode<T>* nodept, std::ostream& os_nodes, bool use=false ) const;   

 private:
  // number of measurements
  unsigned int _nym;
  //! @brief Pointer to measurement data
  const std::list<Data>* _Ym;
  // Chebyshev model of parameters
  std::vector<CV> _CVP;
  // Chebyshev model of outputs
  std::vector<CV> _CVY;
  // Interval bounds of outputs
  std::vector<T> _IY;

  //! @brief Gurobi environment for domain reduction
  GRBEnv* _pGRB;
  //! @brief Chebyshev model environment for output bounding
  CM* _pCM;

  //! @brief vector of operation lists for every output evaluation
  std::vector< std::list<const FFOp*> > _opY;
  //! @brief storage vector for dependents evaluation in interval arithmetic
  std::vector<T> _IEVAL;
  //! @brief storage vector for dependents evaluation in Chebyshev model
  std::vector<CV> _CVEVAL;
  //! @brief vector of independent variables for DAG evaluation
  std::vector<FFVar> _pVAR;
  //! @brief interval bounds of independent variables for DAG evaluation
  std::vector<T> _IVAR;
  //! @brief Chebyshev models of independent variables for DAG evaluation
  std::vector<CV> _CVVAR;


  //! @brief Typedef for cuts
  typedef std::pair<CPPL::dcovector,double> CUT;
  //! @brief Rotate object <a>P</a> into new basis <a>A</a>
  template <typename U>
  std::vector<U> _rotate
    ( const std::vector<U>&P, const CPPL::dgematrix&A ) const;
  //! @brief Rotate object <a>P</a> into new basis <a>A</a>
  template <typename U>
  void _rotate
    ( const std::vector<U>&P, const CPPL::dgematrix&A, std::vector<U>&AxP ) const;
  //! @brief Rotate cutting planes <a>CP</a> into new basis <a>A</a> around point <a>c</a>
  std::vector<CUT> _rotate
    ( const std::vector<CUT>&CP, const CPPL::dgematrix&A ) const;
  //! @brief Perform inclusion test 
  template <typename U>
  typename SetInv<T>::STATUS _inclusion_test
    ( const std::vector<U>&Y ) const;
  //! @brief Perform partial inclusion test (outter/undetermined)
  template <typename U>
  typename SetInv<T>::STATUS _partial_inclusion_test
    ( const std::vector<U>&Y, const std::set<unsigned>&nev) const;   
  //! @brief Construct cutting planes for outer-approximation
  template <typename U>
  std::vector<CUT> _cuttingplanes_outer
    ( const std::vector<U>&Y, const std::vector<T>&P ) const;
  //! @brief Construct cutting planes for outer-approximation
  template <typename U>
  std::vector<CUT> _partial_cuttingplanes_outer
    ( const std::vector<U>&Y, const std::vector<T>&P, const std::set< unsigned > &nev ) const;  
  //! @brief Construct cutting planes for outer-approximation
  std::vector<CUT> _cuttingplanes_parent
    ( const CPPL::dgematrix&A, const std::list<typename SetInvNode<T>::PARENTNODE>&L ) const;
  //! @brief Construct cutting planes for inner-approximation
  template <typename U>
  std::vector<CUT> _cuttingplanes_inner
    ( const std::vector<U>&Y, const std::vector<T>&P ) const;
  //! @brief Perform node outer-reduction using LP
  typename SetInv<T>::STATUS _reduction_outer
    ( const std::vector<CUT>&CP, std::vector<T>&P, double&maxred );
  //! @brief Perform node inner-reduction using constraint propagation
  typename SetInv<T>::STATUS _reduction_inner
    ( const std::vector<CUT>&CP, std::vector<T>&P, double&maxred ) const;
  //! @brief Subproblem assessment using interval arithmetic
  typename SetInv<T>::STATUS _assess_IA
    ( SetInvNode<T>*node );
  //! @brief Subproblem assessment using chebyshev arithmetic
  typename SetInv<T>::STATUS _assess_CM
    ( SetInvNode<T>*node );
  //! @brief Subproblem assessment using chebyshev arithmetic
  typename SetInv<T>::STATUS _assess_CM_loop
    ( std::vector<T>&P, CPPL::dgematrix&A, const bool rotate,
      const std::list<typename SetInvNode<T>::PARENTNODE>&L );
  //! @brief User-function to subproblem assessment
  typename SetInv<T>::STATUS assess
    ( SetInvNode<T>*node );

  //! @brief Private methods to block default compiler methods
  NLEGPE(const NLEGPE&);
  NLEGPE& operator=(const NLEGPE&);
};

template <typename T>
inline void
NLEGPE<T>::output_nodes
( std::ostream&os_open, const bool PROJ, const unsigned DPREC ) const
{
  return SetInv<T>::output_nodes( os_open, PROJ, DPREC );
}

template <typename T>
inline const typename SetInv<T>::t_Nodes&
NLEGPE<T>::open_nodes
() const
{
  return SetInv<T>::open_nodes();
}

template <typename T>
inline void
NLEGPE<T>::output_clusters
( std::ostream&os_clus, const double RTOL, const double ATOL,
  const unsigned DPREC ) const
{
  return SetInv<T>::output_clusters( os_clus, RTOL, ATOL, DPREC );
}

template <typename T>
inline void
NLEGPE<T>::output_node
( const SetInvNode<T>* nodept, std::ostream& os_nodes, bool use ) const
{
  if( use ){
	  for( unsigned i=0; i<_np ; ++i ){
		os_nodes << "I("<<mc::Op<T>::l(nodept->P()[i]) <<","<< mc::Op<T>::u(nodept->P()[i]) << "), ";
	  }
	  os_nodes <<"Node #" << nodept->index() <<"\n"; 
  }
  else{
	  for( unsigned i=0; i<_np ; ++i ){
	  	os_nodes << mc::Op<T>::l(nodept->P()[i]) <<"\t"<< mc::Op<T>::u(nodept->P()[i]) << "\t" ;
	  }
	  os_nodes <<"Node #" << nodept->index() <<"\n"; 
  }
}

template <typename T>
inline std::list<T*>
NLEGPE<T>::clusters
( const double rtol, const double atol, std::ostream&os ) const
{
  return SetInv<T>::clusters( rtol, atol, os );
}

template <typename T>
inline double
NLEGPE<T>::solve
( const T*P, const std::list< Data >&Ym, std::ostream&os )
{
  // Variables
  SetInv<T>::variables( _np, P );
  _CVP.resize(_np);

  // Measurement data
  _Ym  = &Ym;
  _CVY.resize( Ym.size() );
  _IY.resize( Ym.size() );

  // Keep track of execution times
  stats.BOUND = stats.OUTRED = stats.INRED = stats.ROTATE = stats.SETINV = 0.;

  // Initialize internal variables
  _pVAR.resize(_np+1);
  for( unsigned i=0; i<_np; i++ ) _pVAR[i] = _pP[i];
  _pVAR[_np] = *_pT;

  _opY.resize(_ny);
  unsigned nyop = 0;
  for( unsigned i=0; i<_ny; i++ ){
    _opY[i] = _pDAG->subgraph( 1, _pY+i );
    if( nyop <  _opY[i].size() ) nyop = _opY[i].size();
  }

  switch( options.OUTPUT_BOUND ){
  case Options::IA:
    delete _pCM; _pCM = 0;
    _CVEVAL.clear();
    _CVVAR.clear();
    _IEVAL.resize( nyop );
    _IVAR.resize( _pVAR.size() );
    break;

  case Options::CM:
    if( !_pGRB ) _pGRB = new GRBEnv;
    delete _pCM; _pCM  = new CM( _np, options.CM_ORDER );
    _pCM->options = options.CMODEL;
    _CVEVAL.resize( nyop );
    _CVVAR.resize( _pVAR.size() );
    _IEVAL.clear();
    _IVAR.clear();
    break;
  }

  // Run set-inversion algorithm
  stats.SETINV = -time();
  SetInv<T>::options = options.SETINV;
  double res = SetInv<T>::solve( os );
  stats.SETINV += time();
  return res;
}

template <typename T>
inline typename SetInv<T>::STATUS
NLEGPE<T>::assess
( SetInvNode<T>*node )
{ 
  try{ 
    switch( options.OUTPUT_BOUND ){
      // Output bounding using interval arithmetic
      case Options::IA:
        return _assess_IA( node );
      case Options::CM:
        return _assess_CM( node );
      default:
        return SetInv<T>::ABORT;
    }
  }
  catch(...){
#if defined (MC__NLEGPE_DEBUG)
    std::cout << "Failure during assessment of node #" << node->index() << "...\n" << std::endl;
#endif
    return SetInv<T>::FAILURE;
  }
  return SetInv<T>::ABORT; // may not reach this point
}

template <typename T>
inline  typename SetInv<T>::STATUS 
NLEGPE<T>::_assess_IA
( SetInvNode<T>*node )
{
  for( unsigned i=0; i<_np; i++ )
     _IVAR[i] = node->P()[i]; // set current parameter bounds

  typename std::list<Data>::const_iterator ity = _Ym->begin();
  for( unsigned k=0; ity != _Ym->end(); ++ity, k++ ){
    _IVAR[_np] = (*ity).t; // set current measurement point
    _pDAG->eval( _opY[(*ity).iy], _IEVAL.data(), 1, _pY+(*ity).iy, _IY.data()+k,
                 _pVAR.size(), _pVAR.data(), _IVAR.data() );
#if defined (MC__NLEGPE_DEBUG)
  std::cout << "t=" << (*ity).t << "  Y(" << (*ity).iy << ")=" << _IY[k] << std::endl;
#endif
  }

  return _inclusion_test( _IY );
}

/*
TO DO:
1) ADD CUTS FOR PARENT NODES WITH SELECTION OF RELEVANT BOUNDS ONLY
2) ADD ROTATION THRESHOLD WHEN CMODEL HAS CONVERGED
*/

template <typename T>
inline typename SetInv<T>::STATUS 
NLEGPE<T>::_assess_CM
( SetInvNode<T>*node )
{
  // Local copy of bounds and coordinates
  std::vector<T> P = node->P();
  CPPL::dgematrix A = node->A();
#if defined (MC__NLEGPE_DEBUG)
  std::cout << "\nInitial Box and Rotation matrix:\n";
  for( unsigned i=0; i<_np; i++ ) std::cout << P[i] << std::endl;
  std::cout << A << std::endl;
#endif

  // Process node w/o rotation
  typename SetInv<T>::STATUS status = _assess_CM_loop( P, A, false, node->parents() );
  if( status != SetInv<T>::UNDETERMINED ) return status;

  if( options.ROTATION_USE ){
    // Process node w/ rotation
    std::vector<T> P0 = P;
    CPPL::dgematrix A0 = A;
    status = _assess_CM_loop( P, A, true, node->parents() );
    if( status != SetInv<T>::UNDETERMINED ) return status;
    // Update node
    if( SetInv<T>::volume(P) < SetInv<T>::volume(P0) ){
      std::cout << "rotation successful!\n";
      node->parents().push_back( std::make_pair(node->P(),node->A()) );
      node->A() = A;
    }
    else
      P = P0;
  }
  node->P() = P;
  return status;
}

template <typename T>
inline typename SetInv<T>::STATUS 
NLEGPE<T>::_assess_CM_loop
( std::vector<T>&P, CPPL::dgematrix&A, const bool rotate,
  const std::list<typename SetInvNode<T>::PARENTNODE>&L )
{
  typename SetInv<T>::STATUS status;
  double maxredout = options.OUTRED_THRES;
  
  
  for( unsigned iredout=0; ; iredout++ ){
//    bool all_evaluated = true;
//    std::set< unsigned > non_evaluated;
    // Chebyshev model of parameters in rotated coordinates
    for( unsigned i=0; i<_np; i++ ) _CVP[i].set( _pCM, i, P[i] );
    _rotate( _CVP, CPPL::t(A), _CVVAR ); // set current parameter bounds
    // Chebyshev model of outputs in rotated coordinates
    typename std::list<Data>::const_iterator ity = _Ym->begin();
    for( unsigned k=0; ity != _Ym->end(); ++ity, k++ ){
      _CVVAR[_np] = (*ity).t; // set current measurement point
	  try{
	  	_pDAG->eval( _opY[(*ity).iy], _CVEVAL.data(), 1, _pY+(*ity).iy, _CVY.data()+k,
	                 _pVAR.size(), _pVAR.data(), _CVVAR.data() );
	  }
	  catch(...) // Catch exceptions in the function evaluation. 
	  { 
		*(_CVY.data()+k) = T( -1e50 , 1e50  );
//		all_evaluated = false;
//		non_evaluated.insert( (*ity).iy ); // If CVY.data()[ (*ity).iy  ] has not been evaluated, its index is added to non_evaluated.
	  }
#if defined (MC__NLEGPE_DEBUG)
    std::cout << "t=" << (*ity).t << "  Y(" << (*ity).iy << ")=" << _CVY[k] << std::endl;
#endif
    }

//    if(!all_evaluated) status = _partial_inclusion_test( _CVY, non_evaluated ); // perform partial test if at least one element has not been evaluated 
//    else status = _inclusion_test( _CVY );
    
    status = _inclusion_test( _CVY );
    
    if( status != SetInv<T>::UNDETERMINED
     || maxredout < options.OUTRED_THRES
     || iredout >= options.OUTRED_MAX ) break;

    std::vector<CUT> CPout = _cuttingplanes_outer( _CVY, P ) ;
    // Generate cutting planes from Chebyshev model linearization
//    if(all_evaluated) CPout = _cuttingplanes_outer( _CVY, P );
//    else              CPout = _partial_cuttingplanes_outer( _CVY, P, non_evaluated ); // cutting planes only for evaluated CVY's
    
//    std::vector<CUT> CPin  = _cuttingplanes_inner( _CVY, P );

    // Apply rotation
//    if( rotate && !iredout ){
//      // Construct extra rotation matrix Q (only first outer-reduction iteration)
//      CPPL::dgematrix Q = ellintersection_ia( CPout ).eigQ().second;
//#if defined (MC__NLEGPE_DEBUG)
//      std::cout << "\nNew rotation matrix:\n" << Q << std::endl;
//#endif
//      // Add cutting planes of previous nodes
//      std::vector<CUT> CPpar = _cuttingplanes_parent( A, L );
//      CPout.insert( CPout.end(), CPpar.begin(), CPpar.end() );
//      // Rotate variable bounds and cutting planes
//      P = _rotate( P, Q );
//#if defined (MC__NLEGPE_DEBUG)
//      std::cout << "\nNew bounds (w/o centering):\n";
//      for( unsigned i=0; i<_np; i++ ) std::cout << P[i] << std::endl;
//#endif
//      CPout = _rotate( CPout, Q );
//      CPin  = _rotate( CPin,  Q );
//      // Update rotation matrix
//      A = Q*A;
//    }
//    int dum; std::cin >> dum; 

    // Tighten variable bounds
//    bool testflag = false;
//    if(!all_evaluated && Op<T>::l(P[0]) <= 4.379700e-04 &&  Op<T>::u(P[0]) >= 4.379700e-04  ){
//	    testflag = true;
//	    std::cout << " Before domain reduction: \n";
//		for( unsigned i=0; i<_np; i++ ) std::cout << P[i] << "\t";
//		std::cout << std::endl;
//		for( unsigned i=0; i<_ny; i++  ) std::cout << _CVY[i].B() << "\t";
//		std::cout << std::endl;
//		std::cout << "Press [ Enter ] to continue..."; std::cin.clear(); std::cin.ignore(); 
//    }
    status = _reduction_outer( CPout, P, maxredout ); // Nothing to change here!
//    if( testflag ){
//	    std::cout << " After domain reduction: \n";
//		for( unsigned i=0; i<_np; i++ ) std::cout << P[i] << "\t";
//		std::cout << std::endl;
//		for( unsigned i=0; i<_ny; i++  ) std::cout << _CVY[i].B() << "\t";
//		std::cout << std::endl;
//		std::cout << "Press [ Enter ] to continue..."; std::cin.clear(); std::cin.ignore(); 
//    }
    
    if( status != SetInv<T>::UNDETERMINED ) break;
//    double maxredin = options.INRED_THRES;
//    for( unsigned iredin=0; options.INRED_MAX; iredin++ ){
//      status = _reduction_inner( CPin, P, maxredin );
//      //std::cout << maxredin << std::endl << iredin << " inner red.\n";
//      //int dum; std::cin >> dum; 
//      if( status != SetInv<T>::UNDETERMINED 
//       || maxredin < options.INRED_THRES
//       || iredin >= options.INRED_MAX ) break;
//    }
//    if( status != SetInv<T>::UNDETERMINED ) break;
  }     
  return status;
}

template <typename T>
template <typename U>
inline typename SetInv<T>::STATUS 
NLEGPE<T>::_inclusion_test
( const std::vector<U>&Y ) const
{
  // Function enclosure and inclusion tests
  bool flag1 = true, flag2 =true;
  typename std::list<Data>::const_iterator ity = _Ym->begin();
  for( unsigned i=0; flag2 && ity != _Ym->end(); ++ity, i++ ){
    if ( Op<U>::l(Y[i]) > (*ity).yu
      || Op<U>::u(Y[i]) < (*ity).yl )  // no intersection
        flag2 = false;
    else if ( Op<U>::l(Y[i]) < (*ity).yl
           || Op<U>::u(Y[i]) > (*ity).yu )
        flag1 =false;
  }
  if (flag1 && flag2)
    return SetInv<T>::INNER;
  if (!flag2)
    return SetInv<T>::OUTER;
  return SetInv<T>::UNDETERMINED; 
}

template <typename T>
template <typename U>
inline typename SetInv<T>::STATUS 
NLEGPE<T>::_partial_inclusion_test
( const std::vector<U>&Y, const std::set< unsigned > &nev ) const
{
  // Function enclosure and inclusion tests Since not every inclusion has been computed you can only conclude OUTER or UNDETERMINED
  typename std::list<Data>::const_iterator ity = _Ym->begin();
  for( unsigned i=0; ity != _Ym->end(); ++ity, i++ ){
  typename std::set< unsigned >::iterator itk = nev.find( i );  //if i is found in nev then Y[i] doesnt exist. 
    if( itk == nev.end()  ){ // if i is not in nev, nev.find(i) returns nev.end()
		if ( Op<U>::l(Y[i]) > (*ity).yu
		  || Op<U>::u(Y[i]) < (*ity).yl ){  // no intersection
        	return SetInv<T>::OUTER;
 		}
     }    
  }
  return SetInv<T>::UNDETERMINED; 
}

template <typename T>
inline typename SetInv<T>::STATUS 
NLEGPE<T>::_reduction_inner
( const std::vector<CUT>&CP, std::vector<T>&P, double&maxred ) const
{
  typename SetInv<T>::STATUS status = SetInv<T>::UNDETERMINED;
  maxred = 0.;
  const std::vector<T> P0 = P;
#if defined (MC__NLEGPE_DEBUG)
  std::cout << "\n Box before inner domain reduction:\n";
  for( unsigned i=0; i<_np; i++ ) cout << P[i] << endl;
#endif

  // Intersect linear inner relaxation constraints
  for( unsigned i=0; i<_np && status==SetInv<T>::UNDETERMINED; i++ ){
    typename std::vector<CUT>::const_iterator it = CP.begin();
    bool unbnd = false, PinLset = false, PinUset = false; 
    double PinL(0.), PinU(0.);
    for( unsigned k=0; it != CP.end(); ++it, k++ ){
      if( isequal( (*it).first(i), 0. ) ){ unbnd = true; break; }
      T Pbnd = (*it).second;
      for( unsigned j=0; j<_np; j++ )
        if( i != j ) Pbnd -= (*it).first(j)*P[j];
      Pbnd /= (*it).first(i);
      if( (*it).first(i)>0 ){
        PinU = PinUset? std::max( PinU, Op<T>::u(Pbnd) ): Op<T>::u(Pbnd), PinUset = true;
#if defined (MC__NLEGPE_DEBUG)
        std::cout << "Cut #" << k << "  PinU = " << PinU << std::endl;
#endif
      }
      else{
        PinL = PinLset? std::min( PinL, Op<T>::l(Pbnd) ): Op<T>::l(Pbnd), PinLset = true;
#if defined (MC__NLEGPE_DEBUG)
        std::cout << "Cut #" << k << "  PinL = " << PinL << std::endl;
#endif
      }
    }
    // Variable range is unbounded or outer parts intersect - cannot conclude
    if( unbnd || PinL < PinU )
      continue;
    // Outer parts on both sides w/o intersecting - no inner reduction possible
    else if( PinL < Op<T>::u(P[i]) && PinU > Op<T>::l(P[i]) )
      continue;
    // Outer parts on neither sides - no update
    else if( PinL >= Op<T>::u(P[i]) && PinU <= Op<T>::l(P[i]) )
      continue;
    // Outer part on the left-side only
    else if( PinL >= Op<T>::u(P[i]) && PinU < Op<T>::u(P[i]) )
      P[i] = T( Op<T>::l(P[i]), PinU );
    // Outer part on the right-side only
    else if( PinU <= Op<T>::l(P[i]) && PinL > Op<T>::l(P[i]) )
      P[i] = T( PinL, Op<T>::u(P[i]) );
    maxred = std::max( 1.-Op<T>::diam(P[i])/Op<T>::diam(P0[i]), maxred );
  }

#if defined (MC__NLEGPE_DEBUG)
  std::cout << "\n Box after inner domain reduction:\n";
  for( unsigned i=0; i<_np; i++ ) cout << P[i] << endl;
#endif
  return status;
}

template <typename T>
inline typename SetInv<T>::STATUS 
NLEGPE<T>::_reduction_outer
( const std::vector<CUT>&CP, std::vector<T>&P, double&maxred )
{
  maxred = 0.;
  const std::vector<T> P0 = P;
  std::vector<GRBVar> var(_np);

  try{
    // Initialize domain reduction problem
    GRBModel LPred = GRBModel( *_pGRB );
    LPred.getEnv().set( GRB_IntParam_OutputFlag, 0 );  
    LPred.getEnv().set( GRB_DoubleParam_FeasibilityTol, options.OUTRED_TOL );
    LPred.getEnv().set( GRB_DoubleParam_OptimalityTol,  options.OUTRED_TOL );

    // Declare variables
    for( unsigned i=0; i<_np; i++ )
      var[i] = LPred.addVar( Op<T>::l( P[i] ), Op<T>::u( P[i] ), 0., GRB_CONTINUOUS );
    LPred.update();
#if defined (MC__NLEGPE_DEBUG)
    cout<<"\n Box before outer domain reduction:\n";
    for( unsigned i=0; i<_np; i++ ) cout << P[i] << endl;
#endif

    // Add linear outer relaxation constraints
    typename std::vector<CUT>::const_iterator it = CP.begin();
    for( ; it != CP.end(); ++it ){
      GRBLinExpr expr;
      for( unsigned i=0; i<_np; i++ )
        expr += (*it).first(i) * var[i];
      double rhs = (*it).second;        
      LPred.addConstr( expr, GRB_LESS_EQUAL, rhs ); 
    }

    // Solve domain reduction problems
    for( unsigned j=0; j<_np; j++ ){ // _np variables
      GRBLinExpr obj( var[j], 1 );

      for( unsigned i=0; i<2; i++ ){ // 1: max; 0: min
        LPred.setObjective( obj, i? GRB_MAXIMIZE: GRB_MINIMIZE );
        LPred.optimize();
        switch( LPred.get( GRB_IntAttr_Status ) ){

          case GRB_OPTIMAL:
            if( i && var[j].get(GRB_DoubleAttr_X) < Op<T>::u( P[j] ) )
                P[j] = T( Op<T>::l( P[j] ), var[j].get(GRB_DoubleAttr_X) );
            else if( !i && var[j].get(GRB_DoubleAttr_X) > Op<T>::l( P[j] ) )
                P[j] = T( var[j].get(GRB_DoubleAttr_X), Op<T>::u( P[j] ) );
            break;

          case GRB_INFEASIBLE:
          case GRB_INF_OR_UNBD:
#if defined (MC__NLEGPE_DEBUG)
            std::cout << "Infeasible LP!\n";
#endif
            return SetInv<T>::OUTER;

          default:
#if defined (MC__NLEGPE_DEBUG)
            std::cout << "LP solution failed!\n";
#endif
            return SetInv<T>::FAILURE;
        }
      }

      maxred = std::max( 1.-Op<T>::diam(P[j])/Op<T>::diam(P0[j]), maxred );
    }

#if defined (MC__NLEGPE_DEBUG)
    cout<<"\n Box after outer domain reduction:\n";
    for( unsigned i=0; i<_np; i++ ) cout << P[i] << endl;
#endif
    return SetInv<T>::UNDETERMINED;

  } // end: try 

  catch (GRBException e){
    std::cout << "Error code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
  }

  catch (...){
    std::cout << "Exception during optimization" << std::endl;
  }

  return SetInv<T>::FAILURE;
}

template <typename T>
template <typename U>
inline std::vector<typename NLEGPE<T>::CUT>
NLEGPE<T>::_cuttingplanes_outer
( const std::vector<U>&Y, const std::vector<T>&P ) const
{
  std::vector<CUT> CP;
  CPPL::dcovector a(_np);

  // add the two hyperplanes that are generated from each time measurement
  typename std::list<Data>::const_iterator ity = _Ym->begin();
  for( unsigned i=0; ity != _Ym->end(); ++ity, i++ ){
    U Yi = Y[i];
    for( unsigned k=0; k<_np; k++ ) a(k) = Yi.linear( k, true );

    double bU = (*ity).yu - Op<U>::l( Yi );
    for(unsigned k=0;k<_np;k++) bU += a(k) * Op<T>::mid( P[k] );
    CP.push_back( std::make_pair( a, bU ) );

    double bL = -(*ity).yl + Op<U>::u( Yi );  
    for(unsigned k=0;k<_np;k++) bL -= a(k) * Op<T>::mid( P[k] );
    CP.push_back( std::make_pair( -a, bL ) );
  }

  // add the hyperplanes that define the boundaries of the current domain
  for( unsigned i=0; i<_np; i++ ){
    a.zero(); a(i) = 1.;
    CP.push_back( std::make_pair( -a, -Op<T>::l( P[i] ) ) );
    CP.push_back( std::make_pair(  a,  Op<T>::u( P[i] ) ) );
  }

#if defined (MC__NLEGPE_DEBUG)
  std::cout << " Outer Cutting Planes a^T.x <= b in rotated coordinates: " << std::endl;
  std::vector< std::pair<CPPL::dcovector,double> >::const_iterator it = CP.begin();
  for( ; it!=CP.end(); ++it ){
    std::cout << (*it).first;  
    std::cout <<" <=  " << (*it).second << std::endl;
  }
#endif

  return CP;
}


template <typename T>
template <typename U>
inline std::vector<typename NLEGPE<T>::CUT>
NLEGPE<T>::_partial_cuttingplanes_outer
( const std::vector<U>&Y, const std::vector<T>&P, const std::set< unsigned > &nev ) const
{
  std::vector<CUT> CP;
  CPPL::dcovector a(_np);

  // add the two hyperplanes that are generated from each time measurement
  typename std::list<Data>::const_iterator ity = _Ym->begin();
  for( unsigned i=0; ity != _Ym->end(); ++ity, i++ ){
    typename std::set< unsigned >::iterator itk = nev.find( i );     
    if( itk == nev.end() ){ 
    	U Yi = Y[i];
    	for( unsigned k=0; k<_np; k++ ) a(k) = Yi.linear( k, true );

    	double bU = (*ity).yu - Op<U>::l( Yi );
    	for(unsigned k=0;k<_np;k++) bU += a(k) * Op<T>::mid( P[k] );
    	CP.push_back( std::make_pair( a, bU ) );

    	double bL = -(*ity).yl + Op<U>::u( Yi );  
    	for(unsigned k=0;k<_np;k++) bL -= a(k) * Op<T>::mid( P[k] );
    	CP.push_back( std::make_pair( -a, bL ) );
	}
  }

  // add the hyperplanes that define the boundaries of the current domain
  for( unsigned i=0; i<_np; i++ ){
    a.zero(); a(i) = 1.;
    CP.push_back( std::make_pair( -a, -Op<T>::l( P[i] ) ) );
    CP.push_back( std::make_pair(  a,  Op<T>::u( P[i] ) ) );
  }

#if defined (MC__NLEGPE_DEBUG)
  std::cout << " Outer Cutting Planes a^T.x <= b in rotated coordinates: " << std::endl;
  std::vector< std::pair<CPPL::dcovector,double> >::const_iterator it = CP.begin();
  for( ; it!=CP.end(); ++it ){
    std::cout << (*it).first;  
    std::cout <<" <=  " << (*it).second << std::endl;
  }
#endif

  return CP;
}

template <typename T>
inline std::vector<typename NLEGPE<T>::CUT>
NLEGPE<T>::_cuttingplanes_parent
( const CPPL::dgematrix&A, const std::list<typename SetInvNode<T>::PARENTNODE>&L ) const
{
  std::vector<CUT> CP;
#if defined (MC__NLEGPE_DEBUG)
  std::cout << " Parent Cutting Planes a^T.x <= b in rotated coordinates: " << std::endl;
  std::cout << " Current base:\n" << A << std::endl; 
#endif

  // add the hyperplanes that define the boundaries of the parent domains
  typename std::list<typename SetInvNode<T>::PARENTNODE>::const_iterator itl = L.begin();
  for( ; itl != L.end(); ++itl ){
#if defined (MC__NLEGPE_DEBUG)
    std::cout << " Orignal bounds:";
    for( unsigned i=0; i<_np; i++ ) std::cout << "  " << (*itl).first[i];
    std::cout << std::endl;
    std::cout << " Orignal base:\n" << (*itl).second << std::endl; 
#endif
    CPPL::dgematrix AA0T = A*CPPL::t((*itl).second);
    for( unsigned i=0; i<_np; i++ ){
      CP.push_back( std::make_pair( -AA0T.col(i), -Op<T>::l( (*itl).first[i] ) ) );
#if defined (MC__NLEGPE_DEBUG)
      std::cout << -AA0T.col(i);  
      std::cout <<" <=  " << -Op<T>::l( (*itl).first[i] ) << std::endl;
#endif
      CP.push_back( std::make_pair(  AA0T.col(i),  Op<T>::u( (*itl).first[i] ) ) );
#if defined (MC__NLEGPE_DEBUG)
      std::cout << AA0T.col(i);  
      std::cout <<" <=  " << Op<T>::u( (*itl).first[i] ) << std::endl;
#endif
    }
  }

  return CP;
}


template <typename T>
template < typename U >
inline std::vector<typename NLEGPE<T>::CUT>
NLEGPE<T>::_cuttingplanes_inner
( const std::vector<U>&Y, const std::vector<T>&P ) const
{
  std::vector<CUT> CP;
  CPPL::dcovector a(_np);

  // add the two hyperplanes that is generated from each time measurement
  typename std::list<Data>::const_iterator ity = _Ym->begin();
  for( unsigned i=0; ity != _Ym->end(); ++ity, i++ ){
    U Yi = Y[i];
    for( unsigned k=0; k<_np; k++ )
      a(k) = Yi.linear( k, true );

    double bU = (*ity).yl - Op<U>::l( Yi );
    for(unsigned k=0;k<_np;k++) bU += a(k) * Op<T>::mid( P[k] );
    CP.push_back( std::make_pair( a, bU ) );

    double bL = -(*ity).yu + Op<U>::u( Yi );  
    for(unsigned k=0;k<_np;k++) bL -= a(k) * Op<T>::mid( P[k] );
    CP.push_back( std::make_pair( -a, bL ) );
  }

#if defined (MC__NLEGPE_DEBUG)
  std::cout << " Inner Cutting Planes a^T.x <= b in rotated coordinates: " << std::endl;
  std::vector<CUT>::const_iterator it = CP.begin();
  for( ; it!=CP.end(); ++it ){
    std::cout << (*it).first;  
    std::cout <<" <=  " << (*it).second << std::endl;
  }
#endif

  return CP;
}

template <typename T>
template < typename U > inline void
NLEGPE<T>::_rotate
( const std::vector<U>&P, const CPPL::dgematrix&A, std::vector<U>&AxP ) const
{
  for( unsigned i=0; i<_np; i++ ){
    AxP[i] = 0.;
    for( unsigned j=0; j<_np; j++ ){
      AxP[i] += A(i,j)*P[j];
    }
  }
}

template <typename T>
template < typename U > inline std::vector<U>
NLEGPE<T>::_rotate
( const std::vector<U>&P, const CPPL::dgematrix&A ) const
{
  std::vector<U> AxP(_np,0.);
  for( unsigned i=0; i<_np; i++ )
    for( unsigned j=0; j<_np; j++ )
      AxP[i] += A(i,j)*P[j];
  return AxP;
}

template <typename T>
inline std::vector<typename NLEGPE<T>::CUT>
NLEGPE<T>::_rotate
( const std::vector<CUT>&CP, const CPPL::dgematrix&A ) const
{
  std::vector<CUT> CProt;
  CPPL::dgematrix I(_np,_np); I.identity();
  std::vector<CUT>::const_iterator it = CP.begin();
  for( ; it!=CP.end(); ++it ) // aT·y <= b && y = AT·z ==> aT·AT·z <= b
    CProt.push_back( std::make_pair( A*(*it).first, (*it).second ) );

#if defined (MC__NLEGPE_DEBUG)
  std::cout << " Rotated Cutting Planes a^T.x <= b: " << std::endl;
  for( it = CProt.begin(); it!=CProt.end(); ++it ){
    std::cout << (*it).first;  
    std::cout <<" <=  " << (*it).second << std::endl;
  }
#endif
  return CProt;
}


} // namespace mc
#endif

