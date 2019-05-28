#!/usr/bin/gnuplot
export DATA=test8

gnuplot -persist <<-EOFMarker

NX = 9  # Number of states
NQ = 0  # Number of state quadratures

unset key
#set terminal png size 2000,1500 enhanced font "Helvetica,30"
set term post eps enh color solid 6

set output sprintf("%s_STA.eps","$DATA")
set multiplot layout 3,4 title "States"
set style fill transparent solid 0.65 noborder
do for [i = 1:NX]{
  #set xrange [0:1e-4]
  #set yrange [-1:1]
  set xlabel "t"
  set ylabel sprintf("x_%d",i)
  lb = 2*i 
  ub = 2*i+1
  plot sprintf("%s_APPROX_STA.dat","$DATA") u 1:lb:ub w filledcu lt rgb "cyan", \
       sprintf("%s_DINEQI_STA.dat","$DATA") u 1:lb w lines lt rgb "blue", \
       sprintf("%s_DINEQI_STA.dat","$DATA") u 1:ub w lines lt rgb "blue", \
       sprintf("%s_DINEQPM_STA.dat","$DATA") u 1:lb w lines lt rgb "red", \
       sprintf("%s_DINEQPM_STA.dat","$DATA") u 1:ub w lines lt rgb "red", \
       sprintf("%s_EXPANDI_STA.dat","$DATA") u 1:lb w lines lt rgb "blue" dt 3, \
       sprintf("%s_EXPANDI_STA.dat","$DATA") u 1:ub w lines lt rgb "blue" dt 3, \
       sprintf("%s_EXPANDPM_STA.dat","$DATA") u 1:lb w lines lt rgb "red" dt 3, \
       sprintf("%s_EXPANDPM_STA.dat","$DATA") u 1:ub w lines lt rgb "red" dt 3
}
do for [i = 1:NQ]{
  set xlabel "t"
  set ylabel sprintf("q_%d",i)
  lb = 2*NX+2*i 
  ub = 2*NX+2*i+1
  plot sprintf("%s_APPROX_STA.dat","$DATA") u 1:lb:ub w filledcu lt rgb "cyan", \
       sprintf("%s_DINEQI_STA.dat","$DATA") u 1:lb w lines lt rgb "blue", \
       sprintf("%s_DINEQI_STA.dat","$DATA") u 1:ub w lines lt rgb "blue", \
       sprintf("%s_DINEQPM_STA.dat","$DATA") u 1:lb w lines lt rgb "red", \
       sprintf("%s_DINEQPM_STA.dat","$DATA") u 1:ub w lines lt rgb "red", \
       sprintf("%s_EXPANDI_STA.dat","$DATA") u 1:lb w lines lt rgb "blue" dt 3, \
       sprintf("%s_EXPANDI_STA.dat","$DATA") u 1:ub w lines lt rgb "blue" dt 3, \
       sprintf("%s_EXPANDPM_STA.dat","$DATA") u 1:lb w lines lt rgb "red" dt 3, \
       sprintf("%s_EXPANDPM_STA.dat","$DATA") u 1:ub w lines lt rgb "red" dt 3
}
unset multiplot

EOFMarker

ps2eps -B -f -l $(printf "%s_STA.eps" $DATA)
mv $(printf "%s_STA.eps.eps" $DATA) $(printf "%s_STA.eps" $DATA)
gv $(printf "%s_STA.eps" $DATA) &

