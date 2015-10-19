file = 'test6.out'

set size ratio 1
unset key

set xlabel 'D'
set ylabel 'X1'
plot file u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 w boxxy lt 1 lc 1, \
       '' u (($1+$2)/2.):(($3+$4)/2.) w p ps .5

pause -1

set xlabel 'D'
set ylabel 'X2'
plot file u (($1+$2)/2.):(($5+$6)/2.):1:2:5:6 w boxxy lt 1 lc 1, \
       '' u (($1+$2)/2.):(($5+$6)/2.) w p ps .5

pause -1

set xlabel 'D'
set ylabel 'S1'
plot file u (($1+$2)/2.):(($7+$8)/2.):1:2:7:8 w boxxy lt 1 lc 1, \
       '' u (($1+$2)/2.):(($7+$8)/2.) w p ps .5

pause -1

set xlabel 'D'
set ylabel 'S2'
plot file u (($1+$2)/2.):(($9+$10)/2.):1:2:9:10 w boxxy lt 1 lc 1, \
       '' u (($1+$2)/2.):(($9+$10)/2.) w p ps .5

pause -1

set xlabel 'D'
set ylabel 'Z'
plot file u (($1+$2)/2.):(($11+$12)/2.):1:2:11:12 w boxxy lt 1 lc 1, \
       '' u (($1+$2)/2.):(($11+$12)/2.) w p ps .5

pause -1

set xlabel 'D'
set ylabel 'C'
plot file u (($1+$2)/2.):(($13+$14)/2.):1:2:13:14 w boxxy lt 1 lc 1, \
       '' u (($1+$2)/2.):(($13+$14)/2.) w p ps .5

pause -1

file = 'test6b.out'

set xlabel 'D'
set ylabel 'X1'
plot file u 1:2:3 w filledcurves lt 1, \
     '' u 1:2 w l lt 1, '' u 1:3 w l lt 1
#plot file u 1:2:3 w filledcurves lt 2, \
#     '' u 1:2 w l lt 1, '' u 1:3 w l lt 1 lw 1

pause -1

set xlabel 'D'
set ylabel 'X2'
plot file u 1:4:5 w filledcurves lt 1, \
     '' u 1:4 w l lt 1, '' u 1:5 w l lt 1
#plot file u 1:4:5 w filledcurves lt 2, \
#     '' u 1:4 w l lt 1, '' u 1:5 w l lt 1 lw 1

pause -1

set xlabel 'D'
set ylabel 'S1'
plot file u 1:6:7 w filledcurves lt 1, \
     '' u 1:6 w l lt 1, '' u 1:7 w l lt 1

pause -1

set xlabel 'D'
set ylabel 'S2'
plot file u 1:8:9 w filledcurves lt 1, \
     '' u 1:8 w l lt 1, '' u 1:9 w l lt 1

pause -1

set xlabel 'D'
set ylabel 'Z'
plot file u 1:10:11 w filledcurves lt 1, \
     '' u 1:10 w l lt 1, '' u 1:11 w l lt 1

pause -1

set xlabel 'D'
set ylabel 'C'
plot file u 1:12:13 w filledcurves lt 1, \
     '' u 1:12 w l lt 1, '' u 1:13 w l lt 1

set term post eps enh color 18
set out 'test6.eps'
rep
set term x11
!ps2eps -B -f -l test6.eps
!mv test6.eps.eps test6.eps

pause -1

