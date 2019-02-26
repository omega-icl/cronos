file  = 'test0.out'

set term push
set term post eps enh color 10
set out 'test0.eps'

set size ratio 1
#set style fill transparent solid 0.5 noborder
set key spacing 2
#unset key

set xrange [-1.2:1.2]
set yrange [-1.2:1.2]
set xlabel 'p1'
set ylabel 'p2'
plot file u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 tit 'solution set' w boxxy lt 1 lc 1

set term pop
!ps2eps -B -f -l test0.eps
!mv test0.eps.eps test0.eps
!gv test0.eps &


