#!/usr/bin/gnuplot
export DATA=test0

gnuplot -persist <<-EOFMarker

NP = 1  # Number of parameters
NX = 2  # Number of states
NQ = 1  # Number of state quadratures
NF = 2  # Number of state functions

unset key
#set terminal png size 2000,1500 enhanced font "Helvetica,30"
set term post eps enh color solid 10


set output sprintf("%s_STA.eps","$DATA")
set multiplot layout 2,NX+NQ title "States"
set style fill transparent solid 0.65 noborder
do for [i = 1:NX]{
  set xlabel "t"
  set ylabel sprintf("x%d",i)
  plot sprintf("%s_STA.dat","$DATA") u 1:i+1 w l lt rgb "red" lw 2
}
do for [i = 1:NQ]{
  set xlabel "t"
  set ylabel sprintf("q%d",i)
  plot sprintf("%s_STA.dat","$DATA") u 1:NX+i+1 w l lt rgb "red" lw 2
}
unset multiplot

EOFMarker

ps2eps -B -f -l $(printf "%s_STA.eps" $DATA)
mv $(printf "%s_STA.eps.eps" $DATA) $(printf "%s_STA.eps" $DATA)
gv $(printf "%s_STA.eps" $DATA) &

