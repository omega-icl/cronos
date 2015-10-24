
file1 = 'test2_5+1var.out'

set size ratio 1
unset key

set xlabel 'D'
set ylabel 'X1'
plot file1 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 w boxxy lt 1 lc 1#, \
	 #'clusters.out'     u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 w boxxy lt 1 lc 3 fs solid
pause -1

set xlabel 'D'
set ylabel 'X2'
plot file1 u (($1+$2)/2.):(($5+$6)/2.):1:2:5:6 w boxxy lt 1 lc 1#, \
	 #'clusters.out'     u (($1+$2)/2.):(($5+$6)/2.):1:2:5:6 w boxxy lt 1 lc 3 fs solid
pause -1

set xlabel 'D'
set ylabel 'S1'
plot file1 u (($1+$2)/2.):(($7+$8)/2.):1:2:7:8 w boxxy lt 1 lc 1#, \
	 #'clusters.out'     u (($3+$4)/2.):(($5+$6)/2.):3:4:5:6 w boxxy lt 1 lc 3 fs solid
pause -1

set xlabel 'D'
set ylabel 'S2'
plot file1 u (($1+$2)/2.):(($9+$10)/2.):1:2:9:10 w boxxy lt 1 lc 1#, \
	 #'clusters.out'     u (($3+$4)/2.):(($5+$6)/2.):3:4:5:6 w boxxy lt 1 lc 3 fs solid
pause -1

set xlabel 'D'
set ylabel 'Z'
plot file1 u (($1+$2)/2.):(($11+$12)/2.):1:2:11:12 w boxxy lt 1 lc 1#, \
	 #'clusters.out'     u (($3+$4)/2.):(($5+$6)/2.):3:4:5:6 w boxxy lt 1 lc 3 fs solid
pause -1

set xlabel 'D'
set ylabel 'C'
plot file1 u (($1+$2)/2.):(($13+$14)/2.):1:2:13:14 w boxxy lt 1 lc 1#, \
	 #'clusters.out'     u (($3+$4)/2.):(($5+$6)/2.):3:4:5:6 w boxxy lt 1 lc 3 fs solid
pause -1

set term post eps enh color 18
set out 'undetermined.eps'
rep
set term wxt
!ps2eps -B -f -l undetermined.eps
!mv undetermined.eps.eps undetermined.eps

