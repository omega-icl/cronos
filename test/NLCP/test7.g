file1  = 'test7a_95%_Nm4.out'
file2  = 'test7b_95%_Nm4.out'
file3  = 'test7c_95%_Nm4.out'
file4  = 'test7d_95%_Nm4.out'

set term post eps enh color 10
set out 'test7.eps'

set size ratio 1
#set style fill transparent solid 0.5 noborder
set key spacing 2
#unset key

#set xrange [14:32]
#set yrange [0:1]
set xlabel 'p1'
set ylabel 'p2'
#plot file1 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 tit 'enclosure' w boxxy lt 1 lc 1, \
#     file2 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 tit 'exact' w boxxy lt 1 lc 2, \
#     file3 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 tit 'interval bounds' w boxxy lt 1 lc -1
plot file4 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 tit 'LR confidence' w boxxy lt 1 lc 1, \
     file2 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 tit 'exact' w boxxy lt 1 lc 2, \
     file3 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 tit 'interval bounds' w boxxy lt 1 lc -1
#plot file4 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 tit 'LR confidence' w boxxy lt 1 lc 1

set term wxt
!ps2eps -B -f -l test7.eps
!mv test7.eps.eps test7.eps
!gv test7.eps &


