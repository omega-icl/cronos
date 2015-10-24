///////////////////////////////// EXAMPLE 1 ///////////////////////////////////
#if defined( EX1 )

const unsigned NP = 4;                      // # of variables
const I p[NP] = { I(0.,1.), I(0.,1.), I(0.,1.), I(0.,1.) };
//const I p[NP] = { I(2.,60.), I(0.,1.), I(-30.,-1.), I(0.,0.5) };

const unsigned NY = 1;                      // # of output functions
const unsigned NT = 10;                     // # of Measurement points
const double t[NT] = { 0.75, 1.5,  2.25, 3.,     6.,    9.,   13.,   17.,   21.,   25. };
const double y[NT] = { 7.39, 4.09, 1.74, 0.097, -2.57, -2.71, -2.07, -1.44, -0.98, -0.66 };
const unsigned *iy = 0;
const double dy_abs = 0.5, dy_rel = 0.0;

template <class T>
T output_function( const T*p, const T&t, const unsigned i=0 )
{
  //return p[0]*exp(-p[1]*t) + p[2]*exp(-p[3]*t);
  return (58.*p[0]+2.)*exp(-p[1]*t) + (29.*p[2]-30.)*exp(-0.5*p[3]*t);
}

///////////////////////////////// EXAMPLE 2 ///////////////////////////////////
#elif defined( EX2 )

const unsigned NP = 2;                      // # of variables
const I p[NP] = { I(-3.,-0.), I(-3.,3.) };

const unsigned NY = 1;                      // # of output functions
const unsigned NT = 1;                      // # of measurement points
const double t[NT] = { 0. };
const double y[NT] = { 0. };
const unsigned *iy = 0;
const double dy_abs = 0.5, dy_rel = 0.0;

template <class T>
T output_function( const T*p, const T&t, const unsigned i=0 )
{
  return pow(p[0],4)-pow(p[0],2)+4.*pow(p[1],2);
}

///////////////////////////////// EXAMPLE 3 ///////////////////////////////////
#elif defined( EX3 )
      
const unsigned NP = 2;                      // # of variables
//const I p[NP] = { I(-3.,3.), I(-3.,3.) };
const I p[NP] = { I(-1.41667,0.), I(-3.,0.) };
//const I p[NP] = { I(-1.1,-0.9), I(-.2,0.) };
//const I p[NP] = { I(1.176649,1.205612), I(0.7396235,0.7623077) };

const unsigned NY = 1;                      // # of output functions
const unsigned NT = 1;                      // # of Measurement points
const double t[NT] = { 0. };
const double y[NT] = { 1.5 };
const unsigned *iy = 0;
const double dy_abs = 0.5, dy_rel = 0.0;

template <class T>
T output_function( const T*p, const T&t, const unsigned i=0 )
{
  return  sqr(p[0])+sqr(p[1]);
}

///////////////////////////////// EXAMPLE 4 ///////////////////////////////////
#elif defined( EX4 )

const unsigned NP = 2;                      // # of variables
const I p[NP] = { I(-5.,5.), I(-5.,5.) };

const unsigned NY = 1;                      // # of output functions
const unsigned NT = 3;                      // # of Measurement points
const double t[NT] = { -1., 1. , 2. };
const double y[NT] = { -3., 6.5, 9. };
const unsigned *iy = 0;
const double dy_abs = 2.0, dy_rel = 0.0;

template <class T>
T output_function( const T*p, const T&t, const unsigned i=0 )
{
  return sqr(t*p[0])+t*sqr(p[1])+sin(p[0]+t*p[1]);
}

///////////////////////////////// EXAMPLE 5 ///////////////////////////////////
#elif defined( EX5 )
      
const unsigned NP = 3;                      // # of variables
const I p[NP] = { I(-3.,3.), I(-3.,3.), I(-3.,3.) };
//const I p[NP] = { I(-1.41667,0.), I(-3.,0.), I(-2.,-1.) };
//const I p[NP] = { I(-0.6,-0.4), I(-1.,-0.8), I(-1.,-0.8) };

const unsigned NY = 1;                      // # of output functions
const unsigned NT = 1;                      // # of Measurement points
const double t[NT] = { 0. };
const double y[NT] = { 1.5 };
const unsigned *iy = 0;
const double dy_abs = 0.5, dy_rel = 0.0;

template <class T>
T output_function( const T*p, const T&t, const unsigned i=0 )
{
  return  sqr(p[0])+sqr(p[1])+sqr(p[2]);
}

///////////////////////////////// EXAMPLE 8 ///////////////////////////////////
#elif defined( EX8 )

const unsigned NP = 2;                      // # of variables
double e = std::exp(1.);
const I p[NP] = { I(0.25,1.), I(1.5,2*PI) };

const unsigned NY = 1;                      // # of output functions
const unsigned NT = 1;                      // # of measurement points
const double t[NT] = { 0. };
const double y[NT] = { 0. };
const unsigned *iy = 0;
const double dy_abs = 0.0, dy_rel = 0.0;

template <class T>
T output_function( const T*p, const T&t, const unsigned i=0 )
{
  T f1 = 0.5*sin(p[0]*p[1])-0.25*p[1]/PI-0.5*p[0];
  T f2 = (1-0.25/PI)*(exp(2*p[0])-e)+e*p[1]/PI-2*e*p[0];
  return sqr(f1)+sqr(f2);
}

///////////////////////////////// EXAMPLE 8B ///////////////////////////////////
#elif defined( EX8B )

const unsigned NP = 2;                      // # of variables
double e = std::exp(1.);
const I p[NP] = { I(0.25,1.), I(1.5,2*PI) };

const unsigned NY = 2;                      // # of output functions
const unsigned NT = 2;                      // # of measurement points
const double t[NT] = { 0., 0. };
const double y[NT] = { 0., 0. };
const unsigned iy[NT] = { 0, 1 };
const double dy_abs = 1.e-6, dy_rel = 0.0;

template <class T>
T output_function( const T*p, const T&t, const unsigned i=0 )
{
  switch( i ){
  case 0:
    return 0.5*sin(p[0]*p[1])-0.25*p[1]/PI-0.5*p[0];
  case 1: default:
    return (1-0.25/PI)*(exp(2*p[0])-e)+e*p[1]/PI-2*e*p[0];
  };
}

///////////////////////////////// EXAMPLE 10 ///////////////////////////////////
#elif defined( EX10 )

const unsigned NP = 3;                      // # of variables
const I p[NP] = { I(-1.,1.), I(-10.,10.), I(-10.,10.) };

const unsigned NY = 2;                      // # of output functions
const unsigned NT = 2;                      // # of measurement points
const double t[NT] = { 0., 0. };
const double y[NT] = { 0., 0. };
const unsigned iy[NT] = { 0, 1 };
const double dy_abs = 0.0, dy_rel = 0.0;

template <class T>
T output_function( const T*p, const T&t, const unsigned i=0 )
{
  switch( i ){
  case 0:
    return (p[0]+3.)*p[1] + p[0]*p[2] - 3.*p[0];
  case 1: default:
    return p[0]*p[1] + (p[0]+3.)*p[2];
  }
}

#endif
