NX = 2
NF = 2
unset key
#set terminal png size 2000,1500 enhanced font "Helvetica,30"
set term post eps enh color solid 15

SAMPLINGdata  = 'test1_APPROX_STA.dat'
BOUNDINGdata  = 'test1_DINEQI_STA.dat'
BOUNDINGdata2 = 'test1_DINEQPM_STA.dat'

set output "test1_DINEQ_STA.eps"
set multiplot layout NF,NX title "State Bounds"
set style fill transparent solid 0.65 noborder
do for [i = 1:NX]{
  set title sprintf("x%d",i)
  set xlabel "t"
  set ylabel "x"
  lb = 2*i 
  ub = 2*i + 1
  plot SAMPLINGdata u 1:lb:ub w filledcu lt rgb "cyan", \
       BOUNDINGdata u 1:lb w lines lt rgb "green", \
       BOUNDINGdata u 1:ub w lines lt rgb "green", \
       BOUNDINGdata2 u 1:lb w lines lt rgb "red", \
       BOUNDINGdata2 u 1:ub w lines lt rgb "red"
  }
unset multiplot
!ps2eps -B -f -l test1_DINEQ_STA.eps
!mv test1_DINEQ_STA.eps.eps test1_DINEQ_STA.eps 

SAMPLINGdata  = 'test1_APPROX_SEN.dat'
BOUNDINGdata  = 'test1_DINEQI_SEN.dat'
BOUNDINGdata2 = 'test1_DINEQPM_SEN.dat'

set output "test1_DINEQ_SEN.eps"
set multiplot layout NF,NX title "Sensitivity Bounds"
set style fill transparent solid 0.65 noborder
do for [i = 1:NX]{
  set title sprintf("x%d",i)
  set xlabel "t"
  set ylabel "x"
  lb = 2*i 
  ub = 2*i + 1
  plot SAMPLINGdata u 1:lb:ub w filledcu lt rgb "cyan", \
       BOUNDINGdata u 1:lb w lines lt rgb "green", \
       BOUNDINGdata u 1:ub w lines lt rgb "green", \
       BOUNDINGdata2 u 1:lb w lines lt rgb "red", \
       BOUNDINGdata2 u 1:ub w lines lt rgb "red"
  }
unset multiplot
!ps2eps -B -f -l test1_DINEQ_SEN.eps
!mv test1_DINEQ_SEN.eps.eps test1_DINEQ_SEN.eps 

SAMPLINGdata  = 'test1_APPROX_ADJ.dat'
BOUNDINGdata  = 'test1_DINEQI_ADJ.dat'
BOUNDINGdata2 = 'test1_DINEQPM_ADJ.dat'

set output "test1_DINEQ_ADJ.eps"
set multiplot layout NF,NX title "Adjoint Bounds"
set style fill transparent solid 0.65 noborder
do for [i = 0:NF-1]{
  do for [j = 1:NX]{
    set title sprintf("F%d, l%d",i+1,j)
    set xlabel "t"
    set ylabel "l"
    lb = 2*NX*i + 2*j 
    ub = 2*NX*i + 2*j + 1
    plot SAMPLINGdata u 1:lb:ub w filledcu lt rgb "cyan", \
         BOUNDINGdata u 1:lb w lines lt rgb "green", \
         BOUNDINGdata u 1:ub w lines lt rgb "green", \
         BOUNDINGdata2 u 1:lb w lines lt rgb "red", \
         BOUNDINGdata2 u 1:ub w lines lt rgb "red"
  }
}
unset multiplot
!ps2eps -B -f -l test1_DINEQ_ADJ.eps
!mv test1_DINEQ_ADJ.eps.eps test1_DINEQ_ADJ.eps 

