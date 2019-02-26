file1  = 'test7_bnd.out'
file2  = 'test7_inn.out'
file3  = 'test7_clus.out'

set term push
set term post eps enh color 10
set out 'test7.eps'

set size ratio 1
#set style fill transparent solid 0.5 noborder
set key spacing 2
#unset key

set xrange [-1:7]
set yrange [0:5]
set xlabel 'p1'
set ylabel 'p2'
plot file1 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 tit 'boundary set' w boxxy lt 1 lc 1, \
     file2 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 tit 'inner set' w boxxy lt 1 lc 2, \
     file3 u (($1+$2)/2.):(($3+$4)/2.) tit 'FJ point' w p pt 5 lt 1 lc 7

set term pop
!ps2eps -B -f -l test7.eps
!mv test7.eps.eps test7.eps
!gv test7.eps &


