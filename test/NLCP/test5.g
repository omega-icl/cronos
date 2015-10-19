file1 = 'test5b.out'
file = 'test5b_NCO.out'

set size ratio 1
unset key

set xlabel 'alpha'
set ylabel 'beta'
plot file u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 w boxxy lt 1 lc 1, \
     file1 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 w boxxy lt 1 lc 3
#plot file0 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 w boxxy lt 1 lc 1, \
#     '' u (($1+$2)/2.):(($3+$4)/2.) w p ps .5, \
#     file u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 w boxxy lt 1 lc 3, \
#     '' u (($1+$2)/2.):(($3+$4)/2.) w p ps .5

set term post eps enh color 18
set out 'undetermined.eps'
rep
set term x11
!ps2eps -B -f -l undetermined.eps
!mv undetermined.eps.eps undetermined.eps

pause -1

set xlabel 'alpha'
set ylabel 'kappa'
plot file u (($1+$2)/2.):(($5+$6)/2.):1:2:5:6 w boxxy lt 1 lc 1, \
     file1 u (($1+$2)/2.):(($5+$6)/2.):1:2:5:6 w boxxy lt 1 lc 3#,\
     #'' u (($1+$2)/2.):(($5+$6)/2.) w p ps .5

pause -1

set xlabel 'beta'
set ylabel 'kappa'
plot file u (($3+$4)/2.):(($5+$6)/2.):3:4:5:6 w boxxy lt 1 lc 1, \
     file1 u (($3+$4)/2.):(($5+$6)/2.):3:4:5:6 w boxxy lt 1 lc 3#, \
     #'' u (($3+$4)/2.):(($5+$6)/2.) w p ps .5

pause -1

set xlabel 'alpha'
set ylabel 'K'
plot file u (($1+$2)/2.):(($7+$8)/2.):1:2:7:8 w boxxy lt 1 lc 1, \
     '' u (($1+$2)/2.):(($7+$8)/2.) w p ps .5

pause -1

set xlabel 'beta'
set ylabel 'K'
plot file u (($3+$4)/2.):(($7+$8)/2.):3:4:7:8 w boxxy lt 1 lc 1, \
     '' u (($3+$4)/2.):(($7+$8)/2.) w p ps .5

pause -1

set xlabel 'kappa'
set ylabel 'K'
plot file u (($5+$6)/2.):(($7+$8)/2.):5:6:7:8 w boxxy lt 1 lc 1, \
     '' u (($5+$6)/2.):(($7+$8)/2.) w p ps .5

pause -1

