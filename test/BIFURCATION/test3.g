file1 = 'equilibrium.out'
##file2 = 'stable.out'
##file3 = 'unstable.out'
##file4 = 'undetermined.out'
##file5 =  'error.out'

set term post eps enh color 10
set out 'test3.eps'

set grid
unset key
set size 1,1
set origin 0,0
set multiplot layout 3,3 rowsfirst scale .97,1.
#
set xlabel 'a'
set ylabel 'x'
set autoscale y
set yrange [0:2]
plot file1 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 w boxxy lt 1 lc 1 lw 2#, \
#	 file2 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 w boxxy lt 1 lc 3 lw 2, \
#	 file3 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 w boxxy lt 1 lc 5 lw 2, \
#	 file4 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 w boxxy lt 1 lc 6 lw 2, \
#	 file4 u (($1+$2)/2.):(($3+$4)/2.)		   w points     lc 1 pt 1, \
#	 file5 u (($1+$2)/2.):(($3+$4)/2.)         w points     lc 9 pt 2

##
set xlabel 'a'
set ylabel 'y'
set autoscale y
set yrange [0:2]
plot file1 u (($1+$2)/2.):(($5+$6)/2.):1:2:5:6 w boxxy lt 1 lc 1 lw 2#, \
#     file2 u (($1+$2)/2.):(($5+$6)/2.):1:2:5:6 w boxxy lt 1 lc 3 lw 2, \
#	 file3 u (($1+$2)/2.):(($5+$6)/2.):1:2:5:6 w boxxy lt 1 lc 5 lw 2, \
#	 file4 u (($1+$2)/2.):(($5+$6)/2.):1:2:5:6 w boxxy lt 1 lc 7 lw 2, \
#	 file4 u (($1+$2)/2.):(($5+$6)/2.)         w points     lc 1 pt 1, \
#	 file5 u (($1+$2)/2.):(($5+$6)/2.)         w points     lc 9 pt 2
#
set xlabel 'a'
set ylabel 'z'
#set autoscale y
set yrange [0:2]
plot file1 u (($1+$2)/2.):(($7+$8)/2.):1:2:7:8 w boxxy lt 1 lc 1 lw 2#, \
#	 file2 u (($1+$2)/2.):(($7+$8)/2.):1:2:7:8 w boxxy lt 1 lc 3 lw 2, \
#	 file3 u (($1+$2)/2.):(($7+$8)/2.):1:2:7:8 w boxxy lt 1 lc 5 lw 2, \
#	 file4 u (($1+$2)/2.):(($7+$8)/2.):1:2:7:8 w boxxy lt 1 lc 7 lw 2, \
#	 file4 u (($1+$2)/2.):(($7+$8)/2.)         w points     lc 1 pt 1, \
#	 file5 u (($1+$2)/2.):(($7+$8)/2.)         w points     lc 9 pt 2
#
unset multiplot
set term x11

!ps2eps -B -f -l test3.eps
!mv test3.eps.eps test3.eps

