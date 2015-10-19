file = 'test3.out'

set size ratio 1
unset key

set xlabel 'x1'
set ylabel 'x2'
plot file u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 w boxxy lt 1 lc 1, \
       '' u (($1+$2)/2.):(($3+$4)/2.) w p ps .5

pause -1

set xlabel 'x1'
set ylabel 'x3'
plot file u (($1+$2)/2.):(($5+$6)/2.):1:2:5:6 w boxxy lt 1 lc 1, \
       '' u (($1+$2)/2.):(($5+$6)/2.) w p ps .5

pause -1

set xlabel 'x2'
set ylabel 'x3'
plot file u (($3+$4)/2.):(($5+$6)/2.):3:4:5:6 w boxxy lt 1 lc 1, \
       '' u (($3+$4)/2.):(($5+$6)/2.) w p ps .5

pause -1

file = 'test3b.out'

set xlabel 'x1'
set ylabel 'x2'
plot file u 1:2:3 w filledcurves lt 1, \
     '' u 1:2 w l lt 1, '' u 1:3 w l lt 1
#plot file u 1:2:3 w filledcurves lt 2, \
#     '' u 1:2 w l lt 1, '' u 1:3 w l lt 1 lw 1

pause -1

set xlabel 'x1'
set ylabel 'x3'
plot file u 1:4:5 w filledcurves lt 1, \
     '' u 1:4 w l lt 1, '' u 1:5 w l lt 1
#plot file u 1:4:5 w filledcurves lt 2, \
#     '' u 1:4 w l lt 1, '' u 1:5 w l lt 1 lw 1

set term post eps enh color 18
set out 'test3.eps'
rep
set term x11
!ps2eps -B -f -l test3.eps
!mv test3.eps.eps test3.eps

pause -1

