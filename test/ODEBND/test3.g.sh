#!/usr/bin/gnuplot
export DATA=test3

gnuplot -persist <<-EOFMarker

NP = 1  # Number of parameters
NX = 4  # Number of states
NQ = 1  # Number of state quadratures
NF = 1  # Number of state functions

unset key
#set terminal png size 2000,1500 enhanced font "Helvetica,30"
set term post eps enh color solid 10


set output sprintf("%s_STA.eps","$DATA")
set multiplot layout 2,3 title "States"
set style fill transparent solid 0.65 noborder
do for [i = 1:NX]{
  set xlabel "t"
  set ylabel sprintf("x%d",i)
  lb = 2*i 
  ub = 2*i+1
  plot sprintf("%s_APPROX_STA.dat","$DATA") u 1:lb:ub w filledcu lt rgb "cyan", \
       sprintf("%s_DINEQI_STA.dat","$DATA") u 1:lb w lines lt rgb "blue", \
       sprintf("%s_DINEQI_STA.dat","$DATA") u 1:ub w lines lt rgb "blue", \
       sprintf("%s_DINEQPM_STA.dat","$DATA") u 1:lb w lines lt rgb "red", \
       sprintf("%s_DINEQPM_STA.dat","$DATA") u 1:ub w lines lt rgb "red"
}
do for [i = 1:NQ]{
  set xlabel "t"
  set ylabel sprintf("q%d",i)
  lb = 2*NX+2*i 
  ub = 2*NX+2*i+1
  plot sprintf("%s_APPROX_STA.dat","$DATA") u 1:lb:ub w filledcu lt rgb "cyan", \
       sprintf("%s_DINEQI_STA.dat","$DATA") u 1:lb w lines lt rgb "blue", \
       sprintf("%s_DINEQI_STA.dat","$DATA") u 1:ub w lines lt rgb "blue", \
       sprintf("%s_DINEQPM_STA.dat","$DATA") u 1:lb w lines lt rgb "red", \
       sprintf("%s_DINEQPM_STA.dat","$DATA") u 1:ub w lines lt rgb "red"
}
unset multiplot

EOFMarker

ps2eps -B -f -l $(printf "%s_STA.eps" $DATA)
mv $(printf "%s_STA.eps.eps" $DATA) $(printf "%s_STA.eps" $DATA)
gv $(printf "%s_STA.eps" $DATA) &

