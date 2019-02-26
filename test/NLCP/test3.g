file1 = 'test3_bnd.out'

#set size ratio 1
unset key
set style fill transparent solid 0.3 #noborder

#monitorSize=system("xrandr | awk '/\*/{sub(/x/,\",\");print $1; exit}'")
#set macros
#set terminal pngcairo size 1280,720#size @monitorSize
set term post eps enh color 8

################################################################################
set out 'test3.eps'
set multiplot layout 2,2 columnsfirst

set xlabel 'R'
set ylabel 'C1'
plot file1 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 w boxxy lt 1 lc 1 #,\
#     '' u (($1+$2)/2.):(($3+$4)/2.) w p lc 2 pt 9

set ylabel 'C2'
plot file1 u (($1+$2)/2.):(($5+$6)/2.):1:2:5:6 w boxxy lt 1 lc 1 #,\
#     '' u (($1+$2)/2.):(($5+$6)/2.) w p lc 2 pt 9

unset multiplot
#set term wxt
!ps2eps -B -f -l test3.eps
!mv test3.eps.eps test3.eps
!gv test3.eps

