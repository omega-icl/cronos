file1  = 'test11a_90%_Nm32.out'
file2  = 'test11b_90%_Nm32.out'
file3  = 'test11c_90%_Nm32.out'

set term post eps enh color 10
set out 'test11.eps'

set size ratio 1
#set style fill transparent solid 0.5 noborder
set key spacing 2
#unset key

set xrange [14:32]
set yrange [0:1]
set xlabel 'p1'
set ylabel 'p2'
plot file1 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 tit 'enclosure' w boxxy lt 1 lc 1, \
     file2 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 tit 'exact' w boxxy lt 1 lc 2, \
     file3 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 tit 'interval bounds' w boxxy lt 1 lc -1

set term wxt
!ps2eps -B -f -l test11.eps
!mv test11.eps.eps test11.eps
!gv test11.eps &


