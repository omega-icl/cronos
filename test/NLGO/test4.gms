SETS
  N      time stages  / 1*4 /

VARIABLES
  x(N)   decision variables
  z      cost variable;

x.lo(N) = 1;
x.up(N) = 5;

EQUATIONS
  F       Objective function
  G1      Constraint #1
  G2      Constraint #2;



F..  z =E= (x('1')*x('4'))*(x('1')+x('2')+x('3'))+x('3');
G1.. (x('1')*x('4'))*x('2')*x('3') =G= 25;
G2.. power(x('1'),2)+power(x('2'),2)+power(x('3'),2)+power(x('4'),2) =E= 40;

OPTION NLP = ANTIGONE;
OPTION OPTCR = 1e-5;
OPTION OPTCA = 1e-5;
*OPTION DECIMALS = 8;
*$onecho > snopt.opt
*Major feasibility tolerance 1.0e-10
*$offecho

MODEL TEST4 /ALL/;
SOLVE TEST4 USING NLP MINIMIZING z;

