file1 = 'test3.out'
file2 = 'test3b.out'

#set size ratio 1
unset key

#monitorSize=system("xrandr | awk '/\*/{sub(/x/,\",\");print $1; exit}'")
#set macros
#set terminal pngcairo size 1280,720#size @monitorSize
set term post eps enh color 8

################################################################################
set out 'test3.eps'
set multiplot layout 2,2 rowsfirst

set xlabel 'R'
set ylabel 'C1'
plot file1 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 w boxxy lt 1 lc 1

plot file2 u 1:2:3 w filledcurves lt 1, \
     '' u 1:2 w l lt 1, '' u 1:3 w l lt 1       

set ylabel 'C2'
plot file1 u (($1+$2)/2.):(($5+$6)/2.):1:2:5:6 w boxxy lt 1 lc 1

plot file2 u 1:4:5 w filledcurves lt 1, \
     '' u 1:4 w l lt 1, '' u 1:5 w l lt 1       

unset multiplot
set term x11
!ps2eps -B -f -l test3.eps
!mv test3.eps.eps test3.eps
!gv test3.eps

