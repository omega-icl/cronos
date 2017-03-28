file = 'test6.out'

#set size ratio 1
unset key

#monitorSize=system("xrandr | awk '/\*/{sub(/x/,\",\");print $1; exit}'")
#set macros
#set terminal pngcairo size 1280,720#size @monitorSize
set term post eps enh color 8

################################################################################
set out 'test6.eps'
set multiplot layout 2,3 columnsfirst

set xlabel 'D'
set ylabel 'X1'
plot file u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 w boxxy lt 1 lc 1, \
       '' u (($1+$2)/2.):(($3+$4)/2.) w d #p pt 7 ps .1

set ylabel 'X2'
plot file u (($1+$2)/2.):(($5+$6)/2.):1:2:5:6 w boxxy lt 1 lc 1, \
       '' u (($1+$2)/2.):(($5+$6)/2.) w d

set yrange [:5.5]
set ylabel 'S1'
plot file u (($1+$2)/2.):(($7+$8)/2.):1:2:7:8 w boxxy lt 1 lc 1, \
       '' u (($1+$2)/2.):(($7+$8)/2.) w d

set autoscale y
set ylabel 'S2'
plot file u (($1+$2)/2.):(($9+$10)/2.):1:2:9:10 w boxxy lt 1 lc 1, \
       '' u (($1+$2)/2.):(($9+$10)/2.) w d

set ylabel 'Z'
plot file u (($1+$2)/2.):(($11+$12)/2.):1:2:11:12 w boxxy lt 1 lc 1, \
       '' u (($1+$2)/2.):(($11+$12)/2.) w d

set ylabel 'C'
plot file u (($1+$2)/2.):(($13+$14)/2.):1:2:13:14 w boxxy lt 1 lc 1, \
       '' u (($1+$2)/2.):(($13+$14)/2.) w d

unset multiplot
set term wxt
!ps2eps -B -f -l test6.eps
!mv test6.eps.eps test6.eps
!gv test6.eps

