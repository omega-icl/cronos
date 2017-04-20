!export DATA="test2"
data = "test1"

NX = 2
NF = 2
NP = 1
unset key
#set terminal png size 2000,1500 enhanced font "Helvetica,30"
set term post eps enh color solid 10


set output sprintf("%s_STA.eps",data)
set multiplot layout 2,NX title "States"
set style fill transparent solid 0.65 noborder
do for [i = 0:NX-1]{
  set title sprintf("x%d",i)
  set xlabel "t"
  set ylabel "x"
  plot sprintf("%s_STA.dat",data) u 1:i+2 w lp lt rgb "red" lw 2
}
unset multiplot

!ps2eps -B -f -l $(printf "%s_STA.eps" $DATA)
!mv $(printf "%s_STA.eps.eps" $DATA) $(printf "%s_STA.eps" $DATA)
!gv $(printf "%s_STA.eps" $DATA) &


set output sprintf("%s_ASA.eps",data)
set multiplot layout 2,NX title "Adjoint Sensitivities"
set style fill transparent solid 0.65 noborder
do for [i = 0:NF-1]{
  do for [j = 0:NX-1]{
    set title sprintf("f%d, l%d",i,j)
    set xlabel "t"
    set ylabel "l"
    plot sprintf("%s_ASA%d.dat",data,i) u 1:j+2 w lp lt rgb "red" lw 2
  }
}
unset multiplot

!ps2eps -B -f -l $(printf "%s_ASA.eps" $DATA)
!mv $(printf "%s_ASA.eps.eps" $DATA) $(printf "%s_ASA.eps" $DATA)
!gv $(printf "%s_ASA.eps" $DATA) &


set output sprintf("%s_FSA.eps",data)
set multiplot layout 2,NX title "Forward Sensitivities"
set style fill transparent solid 0.65 noborder
do for [j = 0:NX-1]{
  do for [i = 0:NP-1]{
    set title sprintf("dx%d/dp%d",j,i)
    set xlabel "t"
    set ylabel "xp"
    plot sprintf("%s_FSA%d.dat",data,i) u 1:j+2 w lp lt rgb "red" lw 2
  }
}
unset multiplot

!ps2eps -B -f -l $(printf "%s_FSA.eps" $DATA)
!mv $(printf "%s_FSA.eps.eps" $DATA) $(printf "%s_FSA.eps" $DATA)
!gv $(printf "%s_FSA.eps" $DATA) &

