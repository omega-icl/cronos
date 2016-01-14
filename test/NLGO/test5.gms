SETS
  NX     parameters         / 1*3 /
  NT     irradiance stages  / 1*22 /;

VARIABLES
  x(NX)  decision variables
  z      cost variable;

x.lo(NX) = 0;
x.up('1') = 2e-2;
x.up('1') = 1e0;
x.up('1') = 1e0;

PARAMETERS
 Im50(NT) 'Irradiance values for irradiance growth 50 µE/m2·s'
 / 1  10.9249
   2  18.16
   3  36.46
   4  73.12
   5  29.06
   6  69.3541
   7  40
   8  43.54
   9  102.196
   10 94.78
   11 112.995
   12 142.184
   13 207.9
   14 255.617
   15 310.528
   16 387.566
   17 416.857
   18 611.368
   19 640.827
   20 754.596
   21 971.193
   22 1059.39 /
 Pm50(NT) 'Photosynthetic production for irradiance growth 50 µE/m2·s'
 / 1  0.197 
   2  0.412 
   3  0.51955 
   4  0.609 
   5  0.6448 
   6  0.806 
   7  0.788139 
   8  1.074 
   9  1.2006 
   10 1.3617 
   11 1.6483 
   12 2.00669 
   13 2.7235 
   14 2.706 
   15 2.99324 
   16 3.0478 
   17 3.19136 
   18 3.15773 
   19 2.94322 
   20 2.92659 
   21 2.74999 
   22 2.4824 /
 Im1200(NT) 'Irradiance values for irradiance growth 1200 µE/m2·s'
 / 1  21.98 
   2  47.521 
   3  36.41 
   4  50.97 
   5  98.6 
   6  83.83 
   7  160.676 
   8  193.358 
   9  167.544 
   10 222.11 
   11 313.827 
   12 269.622 
   13 419.289 
   14 610.021 
   15 602.581 
   16 624.616 
   17 668.694 
   18 866.909 
   19 998.867 
   20 1186.03 
   21 1310.83
   22 1347.39 /
 Pm1200(NT) 'Photosynthetic production for irradiance growth 1200 µE/m2·s'
 / 1  0.0897199 
   2  0.41226 
   3  0.627 
   4  0.878 
   5  1.0394 
   6  1.2183  
   7  1.703 
   8  2.4369 
   9  2.705 
   10 3.726 
   11 3.781 
   12 4.1386 
   13 5.82 
   14 6.0222 
   15 6.23695 
   16 6.201 
   17 6.11237 
   18 6.00717 
   19 6.33 
   20 6.315 
   21 6.2628 
   22 6.5496 /;

SCALARS
 tau   / 5.5e-3  /
 K     / 3.57e-2 /
 q50   / 8.2e-2  /
 q1200 / 1.8e-2  /;

EQUATIONS
  F       Objective function;
*  G1      Constraint #1
*  G2      Constraint #2;

F..  z =E= SUM(NT, power(x('1')*Im50(NT)/(1.+tau*x('2')*exp(x('3')*log(q50))*Im50(NT)+K*tau*power(x('2')*exp(x('3')*log(q50))*Im50(NT),2))-Pm50(NT),2)) + SUM(NT, power(x('1')*Im1200(NT)/(1.+tau*x('2')*exp(x('3')*log(q1200))*Im1200(NT)+K*tau*power(x('2')*exp(x('3')*log(q1200))*Im1200(NT),2))-Pm1200(NT),2));

*OPTION NLP = BARON;
OPTION NLP = ANTIGONE;
OPTION OPTCR = 1e-3;
OPTION OPTCA = 1e-3;

MODEL TEST5 /ALL/;
SOLVE TEST5 USING NLP MINIMIZING z;

