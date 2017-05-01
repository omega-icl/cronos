#!/usr/bin/gnuplot
export DATA=test4

gnuplot -persist <<-EOFMarker

NP = 3  # Number of parameters
NX = 6  # Number of states
NQ = 0  # Number of state quadratures
NF = 1  # Number of state functions
M  = 4

unset key
#set terminal png size 2000,1500 enhanced font "Helvetica,30"
set term post eps enh color solid 6


set output sprintf("%s_STA.eps","$DATA")
set multiplot layout M,NX+NQ title "States"
set style fill transparent solid 0.65 noborder
do for [i = 1:NX]{
  set xlabel "t"
  set ylabel sprintf("x%d",i)
  lb = 2*i 
  ub = 2*i+1
  plot sprintf("%s_APPROX_STA.dat","$DATA") u 1:lb:ub w filledcu lt rgb "cyan", \
       sprintf("%s_DINEQPM_STA.dat","$DATA") u 1:lb w lines lt rgb "red", \
       sprintf("%s_DINEQPM_STA.dat","$DATA") u 1:ub w lines lt rgb "red", \
       sprintf("%s_DINEQI_STA.dat","$DATA") u 1:lb w lines lt rgb "blue", \
       sprintf("%s_DINEQI_STA.dat","$DATA") u 1:ub w lines lt rgb "blue"
}
do for [i = 1:NQ]{
  set xlabel "t"
  set ylabel sprintf("q%d",i)
  lb = 2*NX+2*i 
  ub = 2*NX+2*i+1
  plot sprintf("%s_APPROX_STA.dat","$DATA") u 1:lb:ub w filledcu lt rgb "cyan", \
       sprintf("%s_DINEQPM_STA.dat","$DATA") u 1:lb w lines lt rgb "red", \
       sprintf("%s_DINEQPM_STA.dat","$DATA") u 1:ub w lines lt rgb "red", \
       sprintf("%s_DINEQI_STA.dat","$DATA") u 1:lb w lines lt rgb "blue", \
       sprintf("%s_DINEQI_STA.dat","$DATA") u 1:ub w lines lt rgb "blue"
}
unset multiplot


set output sprintf("%s_ASA.eps","$DATA")
set multiplot layout M,NX title "Adjoint Sensitivities"
set style fill transparent solid 0.65 noborder
do for [j = 1:NF]{
  do for [i = 1:NX]{
    set xlabel "t"
    set ylabel sprintf("f%d, l%d",j,i)
    lb = 2*i 
    ub = 2*i+1
    plot sprintf("%s_APPROX_ASA%d.dat","$DATA",j-1) u 1:lb:ub w filledcu lt rgb "cyan", \
         sprintf("%s_DINEQPM_ASA%d.dat","$DATA",j-1) u 1:lb w lines lt rgb "red", \
         sprintf("%s_DINEQPM_ASA%d.dat","$DATA",j-1) u 1:ub w lines lt rgb "red", \
         sprintf("%s_DINEQI_ASA%d.dat","$DATA",j-1) u 1:lb w lines lt rgb "blue", \
         sprintf("%s_DINEQI_ASA%d.dat","$DATA",j-1) u 1:ub w lines lt rgb "blue"
  }
  do for [i = 1:NQ]{
    set xlabel "t"
    set ylabel sprintf("f%d, qp%d",j,i)
    lb = 2*NX+2*i 
    ub = 2*NX+2*i+1
    plot sprintf("%s_APPROX_ASA%d.dat","$DATA",j-1) u 1:lb:ub w filledcu lt rgb "cyan", \
         sprintf("%s_DINEQPM_ASA%d.dat","$DATA",j-1) u 1:lb w lines lt rgb "red", \
         sprintf("%s_DINEQPM_ASA%d.dat","$DATA",j-1) u 1:ub w lines lt rgb "red"#, \
#         sprintf("%s_DINEQI_ASA%d.dat","$DATA",j-1) u 1:lb w lines lt rgb "blue", \
#         sprintf("%s_DINEQI_ASA%d.dat","$DATA",j-1) u 1:ub w lines lt rgb "blue"
  }
}
unset multiplot


set output sprintf("%s_FSA.eps","$DATA")
set multiplot layout M,NX+NQ title "Forward Sensitivities"
set style fill transparent solid 0.65 noborder
do for [j = 1:NP]{
  do for [i = 1:NX]{
    set xlabel "t"
    set ylabel sprintf("dx%d/dp%d",i,j)
    lb = 2*i 
    ub = 2*i+1
    plot sprintf("%s_APPROX_FSA%d.dat","$DATA",j-1) u 1:lb:ub w filledcu lt rgb "cyan", \
         sprintf("%s_DINEQPM_FSA%d.dat","$DATA",j-1) u 1:lb w lines lt rgb "red", \
         sprintf("%s_DINEQPM_FSA%d.dat","$DATA",j-1) u 1:ub w lines lt rgb "red"#, \
#         sprintf("%s_DINEQI_FSA%d.dat","$DATA",j-1) u 1:lb w lines lt rgb "blue", \
#         sprintf("%s_DINEQI_FSA%d.dat","$DATA",j-1) u 1:ub w lines lt rgb "blue"
  }
  do for [i = 1:NQ]{
    set xlabel "t"
    set ylabel sprintf("dq%d/dp%d",i,j)
    lb = 2*NX+2*i 
    ub = 2*NX+2*i+1
    plot sprintf("%s_APPROX_FSA%d.dat","$DATA",j-1) u 1:lb:ub w filledcu lt rgb "cyan", \
         sprintf("%s_DINEQPM_FSA%d.dat","$DATA",j-1) u 1:lb w lines lt rgb "red", \
         sprintf("%s_DINEQPM_FSA%d.dat","$DATA",j-1) u 1:ub w lines lt rgb "red"#, \
#         sprintf("%s_DINEQI_FSA%d.dat","$DATA",j-1) u 1:lb w lines lt rgb "blue", \
#         sprintf("%s_DINEQI_FSA%d.dat","$DATA",j-1) u 1:ub w lines lt rgb "blue"
  }
}
unset multiplot

EOFMarker

ps2eps -B -f -l $(printf "%s_STA.eps" $DATA)
mv $(printf "%s_STA.eps.eps" $DATA) $(printf "%s_STA.eps" $DATA)
gv $(printf "%s_STA.eps" $DATA) &

ps2eps -B -f -l $(printf "%s_ASA.eps" $DATA)
mv $(printf "%s_ASA.eps.eps" $DATA) $(printf "%s_ASA.eps" $DATA)
gv $(printf "%s_ASA.eps" $DATA) &

ps2eps -B -f -l $(printf "%s_FSA.eps" $DATA)
mv $(printf "%s_FSA.eps.eps" $DATA) $(printf "%s_FSA.eps" $DATA)
gv $(printf "%s_FSA.eps" $DATA) &

