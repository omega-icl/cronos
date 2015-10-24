splot 'undetermined.out' u (($1+$2)/2.):(($3+$4)/2.):(($5+$6)/2.) w p pt 7, \
	  'clusters.out'	 u (($1+$2)/2.):(($3+$4)/2.):(($5+$6)/2.) w p pt 7
pause -1

set size ratio 1
unset key
set xrange [-0.2:]
set yrange [-0.02:]

set xlabel 'X'
set ylabel 'S'
plot 'undetermined.out' u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 w boxxy lt 1 lc 2, \
	 'clusters.out'     u (($1+$2)/2.):(($3+$4)/2.) w p pt 7 lc 3
pause -1

set xlabel 'X'
set ylabel 'D'
plot 'undetermined.out' u (($1+$2)/2.):(($5+$6)/2.):1:2:5:6 w boxxy lt 1 lc 2, \
	 'clusters.out'     u (($1+$2)/2.):(($5+$6)/2.) w p pt 7 lc 3
pause -1

set xlabel 'S'
set ylabel 'D'
plot 'undetermined.out' u (($3+$4)/2.):(($5+$6)/2.):3:4:5:6 w boxxy lt 1 lc 2, \
	 'clusters.out'     u (($3+$4)/2.):(($5+$6)/2.) w p pt 7 lc 3
pause -1

set term post eps enh color 18
set out 'undetermined.eps'
rep
set term wxt
!ps2eps -B -f -l undetermined.eps
!mv undetermined.eps.eps undetermined.eps

