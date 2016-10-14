file1 = 'test8.out'
file2 = 'test8b.out'
NS = 5
#set size ratio 1
unset key

#monitorSize=system("xrandr | awk '/\*/{sub(/x/,\",\");print $1; exit}'")
#set macros
#set terminal pngcairo size 1280,720#size @monitorSize
set term post eps enh color 6

################################################################################
set out 'test8_methanol.eps'
set multiplot layout 2,NS columnsfirst
do for [i = 1:NS]{
  set title sprintf("x%d",i)
  set xlabel "D"
  set ylabel sprintf("y_{MeOH,%d}",i)
  lb = 2*i+1 
  ub = 2*i+2
  plot file1 u (($1+$2)/2.):((column(lb)+column(ub))/2.):1:2:lb:ub w boxxy lt 1 lc 2#, \
#          '' u (($1+$2)/2.):((column(lb)+column(ub))/2.) w p ps .5
  lb = 2*i 
  ub = 2*i+1
  plot file2 u 1:lb:ub w filledcurves lt 1, \
          '' u 1:lb w l lt 1, '' u 1:ub w l lt 1       
}
unset multiplot

pause -1 "PAUSED"
   
################################################################################
set out 'test8_methyl-butyrate.eps'
set multiplot layout 2,NS columnsfirst
do for [i = 1:NS]{
  set title sprintf("x%d",i)
  set xlabel "D"
  set ylabel sprintf("y_{MeBut,%d}",i)
  lb = 2*NS+2*i+1 
  ub = 2*NS+2*i+2
  plot file1 u (($1+$2)/2.):((column(lb)+column(ub))/2.):1:2:lb:ub w boxxy lt 1 lc 2#, \
#          '' u (($1+$2)/2.):((column(lb)+column(ub))/2.) w p ps .5
  lb = 2*NS+2*i 
  ub = 2*NS+2*i+1
  plot file2 u 1:lb:ub w filledcurves lt 1, \
          '' u 1:lb w l lt 1, '' u 1:ub w l lt 1       
}
unset multiplot

pause -1 "PAUSED"

################################################################################
set out 'test8_toluene.eps'
set multiplot layout 2,NS columnsfirst
do for [i = 1:NS]{
  set title sprintf("x%d",i)
  set xlabel "D"
  set ylabel sprintf("y_{Tol,%d}",i)
  lb = 4*NS+2*i+1 
  ub = 4*NS+2*i+2
  plot file1 u (($1+$2)/2.):((column(lb)+column(ub))/2.):1:2:lb:ub w boxxy lt 1 lc 2#, \
#          '' u (($1+$2)/2.):((column(lb)+column(ub))/2.) w p ps .5
  lb = 4*NS+2*i 
  ub = 4*NS+2*i+1
  plot file2 u 1:lb:ub w filledcurves lt 1, \
          '' u 1:lb w l lt 1, '' u 1:ub w l lt 1       
}
unset multiplot

pause -1 "PAUSED"

################################################################################

