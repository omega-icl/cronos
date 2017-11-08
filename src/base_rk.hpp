// Copyright (C) 2017-... Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef BASE_RK_H
#define BASE_RK_H

#include <iostream>
#include <iomanip>
#include <cmath>

namespace mc
{

//! @brief C++ class for storing/retreiving coefficients of Runge-Kutta methods of various orders
////////////////////////////////////////////////////////////////////////
//! mc::BASE_RK is a C++ class for storing/retreiving coefficients of
//! (explicit) Runge-Kutta methods of various orders
////////////////////////////////////////////////////////////////////////
class BASE_RK
////////////////////////////////////////////////////////////////////////
{
public:

  typedef std::pair<unsigned int,double> t_elem;
  typedef std::pair<unsigned int,t_elem*> t_row;
  //! @brief Default constructor
  BASE_RK
    ( const unsigned int order=0 ):
    _order(0), _stages(0), _tau(0), _alpha(0), _beta(0,0)
    { set(order); }
  //! @brief Constructor
  virtual ~BASE_RK();

  //! @brief Set RK scheme coefficients for order <a>rkord</a>
  bool set
    ( const unsigned int order );

  //! @brief Get RK scheme order
  unsigned int order() const
    { return _order; };
  //! @brief Get RK scheme stages
  unsigned int stages() const
    { return _stages; };
  //! @brief Get RK scheme tau coefficients
  const double* const tau() const
    { return _tau; };
  //! @brief Get RK scheme alpha coefficients
  const t_row* const alpha() const
    { return _alpha; };
  //! @brief Get RK scheme beta coefficients
  const t_row beta() const
    { return _beta; };
  //! @brief Display RK scheme
  void display 
    ( std::ostream&os=std::cout ) const;

protected:

  //! @brief Order of current RK scheme
  unsigned int _order;
  //! @brief Stages in current RK scheme
  unsigned int _stages;
  //! @brief Pointer to current tau coefficients in RK scheme
  double* _tau;
  //! @brief Pointer to current alpha coefficients in RK scheme
  t_row* _alpha;
  //! @brief Pointer to current beta coefficients in RK scheme
  t_row _beta;
};

////////////////////////////////////////////////////////////////////////
inline
BASE_RK::~BASE_RK()
////////////////////////////////////////////////////////////////////////   
{
  if( !_order ) return;
  delete[] _tau;
  for( unsigned int q=0; q<_stages; q++ ) delete[] _alpha[q].second;
  delete[] _alpha;
  delete[] _beta.second;
}

////////////////////////////////////////////////////////////////////////
inline void
BASE_RK::display
( std::ostream&os ) const
////////////////////////////////////////////////////////////////////////   
{
  const unsigned int iprec = 5;
  os << "\nRUNGE-KUTTA SCHEME OF ORDER " << std::left << _order
     << std::endl << std::endl
     << std::scientific << std::setprecision(iprec) << std::right;
  for( unsigned int q=0; q<_stages; q++ ){
    os << std::setw(iprec+7) << _tau[q] << " |";
    unsigned int pos = 0;
    for( unsigned int k=0; k<_alpha[q].first; k++ ){
      for( unsigned int i=pos; i<(_alpha[q].second)[k].first; i++ )
        os << " " << std::setw(iprec+7) << ".";
      os << " " << std::setw(iprec+7) << (_alpha[q].second)[k].second;
      pos = (_alpha[q].second)[k].first+1;
    }
    for( unsigned int i=pos; i<_stages; i++ )
      os << " " << std::setw(iprec+7) << ".";
    os << std::endl;
  }
  os << std::setfill('-') << std::setw((iprec+9)*(1+_stages)-_stages+1) << ""
     << std::setfill(' ') << std::endl
     << std::setw(iprec+9) << "|";
  unsigned int pos = 0;
  for( unsigned int k=0; k<_beta.first; k++ ){
    for( unsigned int i=pos; i<(_beta.second)[k].first; i++ )
      os << " " << std::setw(iprec+7) << ".";
    os << " " << std::setw(iprec+7) << (_beta.second)[k].second;
    pos = (_beta.second)[k].first+1;
  }
  for( unsigned int i=pos; i<_stages; i++ )
    os << " " << std::setw(iprec+7) << ".";
  os << std::endl;
}

////////////////////////////////////////////////////////////////////////
inline bool
BASE_RK::set
( const unsigned int order )
////////////////////////////////////////////////////////////////////////   
{
  if( !order ) return false;

  // Desired RK scheme already defined
  if( order == _order ) return true;
  _order = order;

  // Resize RK scheme
  delete[] _tau;
  for( unsigned int q=0; q<_stages; q++ ) delete[] _alpha[q].second;
  delete[] _alpha;
  delete[] _beta.second;

  switch( _order ){
  case 1: case 2: case 3: case 4: _stages = _order; break;
  case 5: case 6: _stages = _order+1; break;
  //case 7: case 9: _stages = _order+2; break;
  case 8: _stages = _order+3; break;
  default:
    std::cout << "BASE_RK::set *** Selected RK order not available\n";
    return false;
  }
  _tau = new double[_stages];
  _alpha = new t_row[_stages];

  // Define desired RK scheme
  switch( _order ){
  case 1:
    _tau[0] = 0.;

    _alpha[0].first = 0;
    _alpha[0].second = 0;

    _beta.first = 1;
    _beta.second = new t_elem[1];
    (_beta.second)[0] = std::make_pair(0,1.);
    break;

  case 2:
    _tau[0] = 0.;
    _tau[1] = 1.;

    _alpha[0].first = 0;
    _alpha[0].second = 0;
    _alpha[1].first = 1;
    _alpha[1].second = new t_elem[1];
    (_alpha[1].second)[0] = std::make_pair(0,1.);

    _beta.first = 2;
    _beta.second = new t_elem[2];
    (_beta.second)[0] = std::make_pair(0,1./2.);
    (_beta.second)[1] = std::make_pair(1,1./2.);
    break;

  case 3:
    _tau[0] = 0.;
    _tau[1] = 1./2.;
    _tau[2] = 1.;

    _alpha[0].first = 0;
    _alpha[0].second = 0;
    _alpha[1].first = 1;
    _alpha[1].second = new t_elem[1];
    (_alpha[1].second)[0] = std::make_pair(0,1./2.);
    _alpha[2].first = 2;
    _alpha[2].second = new t_elem[2];
    (_alpha[2].second)[0] = std::make_pair(0,-1.);
    (_alpha[2].second)[1] = std::make_pair(1,2.);

    _beta.first = 3;
    _beta.second = new t_elem[3];
    (_beta.second)[0] = std::make_pair(0,1./6.);
    (_beta.second)[1] = std::make_pair(1,2./3.);
    (_beta.second)[2] = std::make_pair(2,1./6.);
    break;

  case 4:
    _tau[0] = 0.;
    _tau[1] = 1./2.;
    _tau[2] = 1./2.;
    _tau[3] = 1.;

    _alpha[0].first = 0;
    _alpha[0].second = 0;
    _alpha[1].first = 1;
    _alpha[1].second = new t_elem[1];
    (_alpha[1].second)[0] = std::make_pair(0,1./2.);
    _alpha[2].first = 1;
    _alpha[2].second = new t_elem[1];
    (_alpha[2].second)[0] = std::make_pair(1,1./2.);
    _alpha[3].first = 1;
    _alpha[3].second = new t_elem[1];
    (_alpha[3].second)[0] = std::make_pair(2,1.);

    _beta.first = 4;
    _beta.second = new t_elem[4];
    (_beta.second)[0] = std::make_pair(0,1./6.);
    (_beta.second)[1] = std::make_pair(1,1./3.);
    (_beta.second)[2] = std::make_pair(2,1./3.);
    (_beta.second)[3] = std::make_pair(3,1./6.);
    break;

  // From the paper: Luther H.A.., ``Further explicit 5th-order Runge-Kutta formulas'',
  //                 SIAM Rev., 8(3):374-380, 1966
  // http://www.jstor.org/stable/2028217
  case 5:
    _tau[0] = 0.;
    _tau[1] = 1.;
    _tau[2] = 1.;
    _tau[3] = 1./4.;
    _tau[4] = 1./2.;
    _tau[5] = 3./4.;

    _alpha[0].first = 0;
    _alpha[0].second = 0;
    _alpha[1].first = 1;
    _alpha[1].second = new t_elem[1];
    (_alpha[1].second)[0] = std::make_pair(0,1.);
    _alpha[2].first = 2;
    _alpha[2].second = new t_elem[2];
    (_alpha[2].second)[0] = std::make_pair(0,1./2.);
    (_alpha[2].second)[1] = std::make_pair(1,1./2.);
    _alpha[3].first = 3;
    _alpha[3].second = new t_elem[3];
    (_alpha[3].second)[0] = std::make_pair(0,14./64.);
    (_alpha[3].second)[1] = std::make_pair(1,5./64.);
    (_alpha[3].second)[2] = std::make_pair(2,-3./64.);
    _alpha[4].first = 4;
    _alpha[4].second = new t_elem[4];
    (_alpha[4].second)[0] = std::make_pair(0,-12./96.);
    (_alpha[4].second)[1] = std::make_pair(1,-12./96.);
    (_alpha[4].second)[2] = std::make_pair(2,8./96.);
    (_alpha[4].second)[3] = std::make_pair(3,64./96.);
    _alpha[5].first = 4;
    _alpha[5].second = new t_elem[4];
    (_alpha[5].second)[0] = std::make_pair(1,-9./64.);
    (_alpha[5].second)[1] = std::make_pair(2,5./64.);
    (_alpha[5].second)[2] = std::make_pair(3,16./64.);
    (_alpha[5].second)[3] = std::make_pair(4,36./64.);

    _beta.first = 5;
    _beta.second = new t_elem[5];
    (_beta.second)[0] = std::make_pair(0,7./90.);
    (_beta.second)[1] = std::make_pair(2,7./90.);
    (_beta.second)[2] = std::make_pair(3,32./90.);
    (_beta.second)[3] = std::make_pair(4,12./90.);
    (_beta.second)[4] = std::make_pair(5,32./90.);
    break;
// 
//   // From the paper: Butcher J.C.,  ``On Runge-Kutta processes of high order'',
//   //                 J. Austral. Math. Soc., 4:179-194, 1964
//   // http://journals.cambridge.org/download.php?file=%2FJAZ%2FJAZ4_02%2FS1446788700023387a.pdf
//   case 6:{
//     const double mu = 1./3.;
//     const double lambda = 4.;
// 
//     _tau[0] = 0.;
//     _tau[1] = mu;
//     _tau[2] = 2./3.;
//     _tau[3] = 1./3.;
//     _tau[4] = 1./2.;
//     _tau[5] = 1./2.;
//     _tau[6] = 1.;
// 
//     _alpha[0].first = 0;
//     _alpha[0].second = 0;
//     _alpha[1].first = 1;
//     _alpha[1].second = new t_elem[1];
//     (_alpha[1].second)[0] = std::make_pair(0,mu);
//     _alpha[2].first = 2;
//     _alpha[2].second = new t_elem[2];
//     (_alpha[2].second)[0] = std::make_pair(0,2./3.-2./(9.*mu));
//     (_alpha[2].second)[1] = std::make_pair(1,2./(9.*mu));
//     _alpha[3].first = 3;
//     _alpha[3].second = new t_elem[3];
//     (_alpha[3].second)[0] = std::make_pair(0,5./12.-1./(9.*mu));
//     (_alpha[3].second)[1] = std::make_pair(1,1./(9.*mu));
//     (_alpha[3].second)[2] = std::make_pair(2,-1./12.);
//     _alpha[4].first = 4;
//     _alpha[4].second = new t_elem[4];
//     (_alpha[4].second)[0] = std::make_pair(0,17./16.-3./(8.*mu));
//     (_alpha[4].second)[1] = std::make_pair(1,3./(8.*mu));
//     (_alpha[4].second)[2] = std::make_pair(2,-3./16.);
//     (_alpha[4].second)[3] = std::make_pair(3,-3./8.);
//     _alpha[5].first = 5;
//     _alpha[5].second = new t_elem[5];
//     (_alpha[5].second)[0] = std::make_pair(0,17./16.-3./(8.*mu)+1./(4.*lambda));
//     (_alpha[5].second)[1] = std::make_pair(1,3./(8.*mu));
//     (_alpha[5].second)[2] = std::make_pair(2,-3./16.-3./(4.*lambda));
//     (_alpha[5].second)[3] = std::make_pair(3,-3./8.-3./(2.*lambda));
//     (_alpha[5].second)[4] = std::make_pair(4,2./lambda);
//     _alpha[6].first = 6;
//     _alpha[6].second = new t_elem[6];
//     (_alpha[6].second)[0] = std::make_pair(0,-27./44.+3./(11.*mu));
//     (_alpha[6].second)[1] = std::make_pair(1,-3./(11.*mu));
//     (_alpha[6].second)[2] = std::make_pair(2,63./44.);
//     (_alpha[6].second)[3] = std::make_pair(3,18./11.);
//     (_alpha[6].second)[4] = std::make_pair(4,(4.*lambda-16.)/11.);
//     (_alpha[6].second)[5] = std::make_pair(5,-4.*lambda/11.);
// 
//     _beta.first = 6;
//     _beta.second = new t_elem[6];
//     (_beta.second)[0] = std::make_pair(0,11./120.);
//     (_beta.second)[1] = std::make_pair(2,27./40.);
//     (_beta.second)[2] = std::make_pair(3,27./40.);
//     (_beta.second)[3] = std::make_pair(4,(lambda-8.)/15.);
//     (_beta.second)[4] = std::make_pair(5,-lambda/15.);
//     (_beta.second)[5] = std::make_pair(6,11./120.);
//     break;
//   }

  // From the paper: Butcher J.C.,  ``On Runge-Kutta processes of high order'',
  //                 J. Austral. Math. Soc., 4:179-194, 1964
  // http://journals.cambridge.org/production/action/cjoGetFulltext?fulltextid%3D4910152&ei=4g7QTsTND4Wy8QPrpOHlDw&usg=AFQjCNF2BwP1ZBNRgwiev9lydOvmcVQc7A&sig2=NAnVVAJvIMA2uLYkm5hGxA&cad=rja

  case 6:{
    _tau[0] = 0.;
    _tau[1] = 1./3.;
    _tau[2] = 2./3.;
    _tau[3] = 1./3.;
    _tau[4] = 1./2.;
    _tau[5] = 1./2.;
    _tau[6] = 1.;

    _alpha[0].first = 0;
    _alpha[0].second = 0;
    _alpha[1].first = 1;
    _alpha[1].second = new t_elem[1];
    (_alpha[1].second)[0] = std::make_pair(0,1./3.);
    _alpha[2].first = 1;
    _alpha[2].second = new t_elem[1];
    (_alpha[2].second)[0] = std::make_pair(1,2./3.);
    _alpha[3].first = 3;
    _alpha[3].second = new t_elem[3];
    (_alpha[3].second)[0] = std::make_pair(0,1./12.);
    (_alpha[3].second)[1] = std::make_pair(1,1./3.);
    (_alpha[3].second)[2] = std::make_pair(2,-1./12.);
    _alpha[4].first = 4;
    _alpha[4].second = new t_elem[4];
    (_alpha[4].second)[0] = std::make_pair(0,-1./16.);
    (_alpha[4].second)[1] = std::make_pair(1,9./8.);
    (_alpha[4].second)[2] = std::make_pair(2,-3./16.);
    (_alpha[4].second)[3] = std::make_pair(3,-3./8.);
    _alpha[5].first = 4;
    _alpha[5].second = new t_elem[4];
    (_alpha[5].second)[0] = std::make_pair(1,9./8.);
    (_alpha[5].second)[1] = std::make_pair(2,-3./8.);
    (_alpha[5].second)[2] = std::make_pair(3,-3./4.);
    (_alpha[5].second)[3] = std::make_pair(4,1./2.);
    _alpha[6].first = 5;
    _alpha[6].second = new t_elem[5];
    (_alpha[6].second)[0] = std::make_pair(0,9./44.);
    (_alpha[6].second)[1] = std::make_pair(1,-9./11.);
    (_alpha[6].second)[2] = std::make_pair(2,63./44.);
    (_alpha[6].second)[3] = std::make_pair(3,18./11.);
    (_alpha[6].second)[4] = std::make_pair(5,-16./11.);

    _beta.first = 6;
    _beta.second = new t_elem[6];
    (_beta.second)[0] = std::make_pair(0,11./120.);
    (_beta.second)[1] = std::make_pair(2,27./40.);
    (_beta.second)[2] = std::make_pair(3,27./40.);
    (_beta.second)[3] = std::make_pair(4,-4./15.);
    (_beta.second)[4] = std::make_pair(5,-4./15.);
    (_beta.second)[5] = std::make_pair(6,11./120.);
    break;
  }

  // From the paper: Cooper G.J., Verner J.H., ``Some Explicit Runge-Kutta Methods of High Order'',
  //                 SIAM J. Numer. Anal., 9(3):389-405, 1972
  // http://www.jstor.org/stable/2156139
  case 8:
    _tau[0] = 0.;
    _tau[1] = 1./2.;
    _tau[2] = 1./2.;
    _tau[3] = (7.-std::sqrt(21.))/14.;
    _tau[4] = (7.-std::sqrt(21.))/14.;
    _tau[5] = 1./2.;
    _tau[6] = (7.+std::sqrt(21.))/14.;
    _tau[7] = (7.+std::sqrt(21.))/14.;
    _tau[8] = 1./2.;
    _tau[9] = (7.-std::sqrt(21.))/14.;
    _tau[10] = 1.;

    _alpha[0].first = 0;
    _alpha[0].second = 0;
    _alpha[1].first = 1;
    _alpha[1].second = new t_elem[1];
    (_alpha[1].second)[0] = std::make_pair(0,1./2.);
    _alpha[2].first = 2;
    _alpha[2].second = new t_elem[2];
    (_alpha[2].second)[0] = std::make_pair(0,1./4.);
    (_alpha[2].second)[1] = std::make_pair(1,1./4.);
    _alpha[3].first = 3;
    _alpha[3].second = new t_elem[3];
    (_alpha[3].second)[0] = std::make_pair(0,1./7.);
    (_alpha[3].second)[1] = std::make_pair(1,(-7.+3.*std::sqrt(21.))/98.);
    (_alpha[3].second)[2] = std::make_pair(2,(21.-5.*std::sqrt(21.))/49.);
    _alpha[4].first = 3;
    _alpha[4].second = new t_elem[3];
    (_alpha[4].second)[0] = std::make_pair(0,(11.-std::sqrt(21.))/84.);
    (_alpha[4].second)[1] = std::make_pair(2,(18.-4.*std::sqrt(21.))/63.);
    (_alpha[4].second)[2] = std::make_pair(3,(21.+std::sqrt(21.))/252.);
    _alpha[5].first = 4;
    _alpha[5].second = new t_elem[4];
    (_alpha[5].second)[0] = std::make_pair(0,(5.-std::sqrt(21.))/48.);
    (_alpha[5].second)[1] = std::make_pair(2,(9.-std::sqrt(21.))/36.);
    (_alpha[5].second)[2] = std::make_pair(3,(-231.-14.*std::sqrt(21.))/360.);
    (_alpha[5].second)[3] = std::make_pair(4,(63.+7.*std::sqrt(21.))/80.);
    _alpha[6].first = 5;
    _alpha[6].second = new t_elem[5];
    (_alpha[6].second)[0] = std::make_pair(0,(10.+std::sqrt(21.))/42.);
    (_alpha[6].second)[1] = std::make_pair(2,(-432.-92.*std::sqrt(21.))/315.);
    (_alpha[6].second)[2] = std::make_pair(3,(633.+145.*std::sqrt(21.))/90.);
    (_alpha[6].second)[3] = std::make_pair(4,(-504.-115.*std::sqrt(21.))/70.);
    (_alpha[6].second)[4] = std::make_pair(5,(63.+13.*std::sqrt(21.))/35.);
    _alpha[7].first = 4;
    _alpha[7].second = new t_elem[4];
    (_alpha[7].second)[0] = std::make_pair(0,1./14.);
    (_alpha[7].second)[1] = std::make_pair(4,(14.+3.*std::sqrt(21.))/126.);
    (_alpha[7].second)[2] = std::make_pair(5,(13.+3.*std::sqrt(21.))/63.);
    (_alpha[7].second)[3] = std::make_pair(6,1./9.);
    _alpha[8].first = 5;
    _alpha[8].second = new t_elem[5];
    (_alpha[8].second)[0] = std::make_pair(0,1./32.);
    (_alpha[8].second)[1] = std::make_pair(4,(91.+21.*std::sqrt(21.))/576.);
    (_alpha[8].second)[2] = std::make_pair(5,11./72.);
    (_alpha[8].second)[3] = std::make_pair(6,(-385.+75.*std::sqrt(21.))/1152.);
    (_alpha[8].second)[4] = std::make_pair(7,(63.-13.*std::sqrt(21.))/128.);
    _alpha[9].first = 6;
    _alpha[9].second = new t_elem[6];
    (_alpha[9].second)[0] = std::make_pair(0,1./14.);
    (_alpha[9].second)[1] = std::make_pair(4,1./9.);
    (_alpha[9].second)[2] = std::make_pair(5,(-733.+147.*std::sqrt(21.))/2205.);
    (_alpha[9].second)[3] = std::make_pair(6,(515.-111.*std::sqrt(21.))/504.);
    (_alpha[9].second)[4] = std::make_pair(7,(-51.+11.*std::sqrt(21.))/56.);
    (_alpha[9].second)[5] = std::make_pair(8,(132.-28.*std::sqrt(21.))/245.);
    _alpha[10].first = 6;
    _alpha[10].second = new t_elem[6];
    (_alpha[10].second)[0] = std::make_pair(4,(-42.-7.*std::sqrt(21.))/18.);
    (_alpha[10].second)[1] = std::make_pair(5,(-18.-28.*std::sqrt(21.))/45.);
    (_alpha[10].second)[2] = std::make_pair(6,(-273.+53.*std::sqrt(21.))/72.);
    (_alpha[10].second)[3] = std::make_pair(7,(301.-53.*std::sqrt(21.))/72.);
    (_alpha[10].second)[4] = std::make_pair(8,(28.+28.*std::sqrt(21.))/45.);
    (_alpha[10].second)[5] = std::make_pair(9,(49.+7.*std::sqrt(21.))/18.);

    _beta.first = 5;
    _beta.second = new t_elem[5];
    (_beta.second)[0] = std::make_pair(0,1./20.);
    (_beta.second)[1] = std::make_pair(7,49./180.);
    (_beta.second)[2] = std::make_pair(8,16./45.);
    (_beta.second)[3] = std::make_pair(9,49./180.);
    (_beta.second)[4] = std::make_pair(10,1./20.);
    break;

  default: return false;
  }

  return true;
}

} // namespace mc

#endif
