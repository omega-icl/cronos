NX = 3
NF = 2
unset key
#set terminal png size 2000,1500 enhanced font "Helvetica,30"
se term post eps enh color solid 10

SAMPLINGdata = 'test2_APPROX_STA.dat'
BOUNDIdata   = 'test2_DINEQI_STA.dat'
BOUNDPMdata  = 'test2_DINEQPM_STA.dat'
BOUNDVALdata = 'test2_VAL_STA.dat'

set output "test2_STA.eps"
set multiplot layout NF,NX title "State Bounds"
set style fill transparent solid 0.65 noborder
do for [i = 1:NX]{
  set title sprintf("x%d",i)
  set xlabel "t"
  set ylabel "x"
  lb = 2*i 
  ub = 2*i + 1
  plot SAMPLINGdata u 1:lb:ub w filledcu lt rgb "cyan", \
       BOUNDIdata u 1:lb w lines lt rgb "magenta", \
       BOUNDIdata u 1:ub w lines lt rgb "magenta", \
       BOUNDPMdata u 1:lb w lines lt rgb "blue", \
       BOUNDPMdata u 1:ub w lines lt rgb "blue", \
       BOUNDVALdata u 1:lb w lines lt rgb "red", \
       BOUNDVALdata u 1:ub w lines lt rgb "red"
  }
unset multiplot

!ps2eps -B -f -l test2_STA.eps
!mv test2_STA.eps.eps test2_STA.eps
!gv test2_STA.eps

pause 0
