file  = 'test0.out'
file2 = 'test0b.out'
file3 = 'test0b2.out'
file4 = 'test0c.out'

set term post eps enh color 10
set out 'test0.eps'

#set size ratio 1
set style fill transparent solid 0.5 noborder
set key spacing 2
#unset key

set xlabel 'p1'
set ylabel 'p2'
plot file4 u 2:3 tit 'SAMPLING' w p ps .5 pt 6 lc 4, \
     file2 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 tit 'U MAX-MIN' w boxxy lt 1 lc 1, \
     file3 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 tit 'U SIP' w boxxy lt 1 lc 3, \
     file u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 tit 'INNER' w boxxy lt 1 lc 2

set term x11
!ps2eps -B -f -l test0.eps
!mv test0.eps.eps test0.eps

