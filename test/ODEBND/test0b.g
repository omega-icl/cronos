NX = 2
NF = 2
unset key
#set terminal png size 2000,1500 enhanced font "Helvetica,30"
se term post eps enh color solid 10

SAMPLINGdata = 'test0_APPROX_STA.dat'
BOUNDIdata   = 'test0_EXPANDI_STA.dat'
BOUNDPMdata  = 'test0_EXPANDPM_STA.dat'

set output "test0_EXPAND_STA.eps"
set multiplot layout NF,NX title "State Bounds"
set style fill transparent solid 0.65 noborder
do for [i = 1:NX]{
  set title sprintf("x%d",i)
  set xlabel "t"
  set ylabel "x"
  set yrange [0:2]
  lb = 2*i 
  ub = 2*i + 1
  plot SAMPLINGdata u 1:lb:ub w filledcu lt rgb "cyan", \
       BOUNDIdata u 1:lb w lines lt rgb "magenta", \
       BOUNDIdata u 1:ub w lines lt rgb "magenta", \
       BOUNDPMdata u 1:lb w lines lt rgb "blue", \
       BOUNDPMdata u 1:ub w lines lt rgb "blue"
  }
unset multiplot

!ps2eps -B -f -l test0_EXPA_STA.eps
!mv test0_EXPAND_STA.eps.eps test0_EXPAND_STA.eps
!gv test0_EXPAND_STA.eps

