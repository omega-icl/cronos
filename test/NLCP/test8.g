file1 = 'test8.out'
NS = 3
#set size ratio 1
unset key

#monitorSize=system("xrandr | awk '/\*/{sub(/x/,\",\");print $1; exit}'")
#set macros
#set terminal pngcairo size 1280,720#size @monitorSize
set term post eps enh color 6

################################################################################
set out 'test8.eps'
set multiplot layout 3,NS rowsfirst
do for [i = 1:NS]{
  set title sprintf("x%d",i)
  set xlabel "D"
  set ylabel sprintf("y_{MeOH,%d}",i)
  lb = 2*i+1 
  ub = 2*i+2
  plot file1 u (($1+$2)/2.):((column(lb)+column(ub))/2.):1:2:lb:ub w boxxy lt 1 lc 2#, \
#          '' u (($1+$2)/2.):((column(lb)+column(ub))/2.) w p ps .5
}
do for [i = 1:NS]{
  set title sprintf("x%d",i)
  set xlabel "D"
  set ylabel sprintf("y_{MeBut,%d}",i)
  lb = 2*NS+2*i+1 
  ub = 2*NS+2*i+2
  plot file1 u (($1+$2)/2.):((column(lb)+column(ub))/2.):1:2:lb:ub w boxxy lt 1 lc 2#, \
#          '' u (($1+$2)/2.):((column(lb)+column(ub))/2.) w p ps .5
}
do for [i = 1:NS]{
  set title sprintf("x%d",i)
  set xlabel "D"
  set ylabel sprintf("y_{Tol,%d}",i)
  lb = 4*NS+2*i+1 
  ub = 4*NS+2*i+2
  plot file1 u (($1+$2)/2.):((column(lb)+column(ub))/2.):1:2:lb:ub w boxxy lt 1 lc 2#, \
#          '' u (($1+$2)/2.):((column(lb)+column(ub))/2.) w p ps .5
}
unset multiplot
################################################################################

set term wxt
!ps2eps -B -f -l test8.eps
!mv test8.eps.eps test8.eps
!gv test8.eps

