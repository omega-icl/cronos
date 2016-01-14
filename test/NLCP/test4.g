file  = 'test4.out'
file2 = 'test4b.out'

set term post eps enh color 10
set out 'test4.eps'

#set size ratio 1
set style fill transparent solid 0.5 noborder
unset key

set multiplot layout 2,3 columnsfirst margins 0.1,0.9,0.1,0.9 spacing 0.1

set xlabel 'x1'
set ylabel 'x2'
plot file2 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 w boxxy lt 1 lc 1, \
     file u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 w boxxy lt 1 lc 2
 
#pause -1

set xlabel 'x1'
set ylabel 'x3'
plot file2 u (($1+$2)/2.):(($5+$6)/2.):1:2:5:6 w boxxy lt 1 lc 1, \
     file u (($1+$2)/2.):(($5+$6)/2.):1:2:5:6 w boxxy lt 1 lc 2

#pause -1

set xlabel 'x2'
set ylabel 'x3'
plot file2 u (($3+$4)/2.):(($5+$6)/2.):3:4:5:6 w boxxy lt 1 lc 1, \
     file u (($3+$4)/2.):(($5+$6)/2.):3:4:5:6 w boxxy lt 1 lc 2

#pause -1

set xlabel 'x1'
set ylabel 'x4'
plot file2 u (($1+$2)/2.):(($7+$8)/2.):1:2:7:8 w boxxy lt 1 lc 1, \
     file u (($1+$2)/2.):(($7+$8)/2.):1:2:7:8 w boxxy lt 1 lc 2

#pause -1

set xlabel 'x2'
set ylabel 'x4'
plot file2 u (($3+$4)/2.):(($7+$8)/2.):3:4:7:8 w boxxy lt 1 lc 1, \
     file u (($3+$4)/2.):(($7+$8)/2.):3:4:7:8 w boxxy lt 1 lc 2

#pause -1

set xlabel 'x3'
set ylabel 'x4'
plot file2 u (($5+$6)/2.):(($7+$8)/2.):5:6:7:8 w boxxy lt 1 lc 1, \
     file u (($5+$6)/2.):(($7+$8)/2.):5:6:7:8 w boxxy lt 1 lc 2

#pause -1
unset multiplot

set term x11
!ps2eps -B -f -l test4.eps
!mv test4.eps.eps test4.eps

