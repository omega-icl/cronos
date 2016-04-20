NX = 2
NF = 2
NP = 1
unset key
#set terminal png size 2000,1500 enhanced font "Helvetica,30"
se term post eps enh color solid 10

data = 'test1_STA.dat'

set output "test1_STA.eps"
set multiplot layout 2,NX title "States"
set style fill transparent solid 0.65 noborder
do for [i = 1:NX]{
  set title sprintf("x%d",i)
  set xlabel "t"
  set ylabel "x"
  plot data u 1:i+1 w l lt rgb "red" lw 2
}
unset multiplot

!ps2eps -B -f -l test1_STA.eps
!mv test1_STA.eps.eps test1_STA.eps
!gv test1_STA.eps

data  = 'test1_SEN.dat'

set output "test1_SEN.eps"
set multiplot layout 2,NX title "State Sensitivities"
set style fill transparent solid 0.65 noborder
do for [j = 1:NX]{
  do for [i = 0:NP-1]{
    set title sprintf("dx%d/dp%d",j,i+1)
    set xlabel "t"
    set ylabel "xp"
    plot data u 1:NX*i+j+1 w l lt rgb "red" lw 2
  }
}
unset multiplot
!ps2eps -B -f -l test1_SEN.eps
!mv test1_SEN.eps.eps test1_SEN.eps 
!gv test1_SEN.eps

data  = 'test1_ADJ.dat'

set output "test1_ADJ.eps"
set multiplot layout 2,NX title "Adjoints"
set style fill transparent solid 0.65 noborder
do for [i = 0:NF-1]{
  do for [j = 1:NX]{
    set title sprintf("f%d, l%d",i+1,j)
    set xlabel "t"
    set ylabel "l"
    plot data u 1:NX*i+j+1 w l lt rgb "red" lw 2
  }
}
unset multiplot
!ps2eps -B -f -l test1_ADJ.eps
!mv test1_ADJ.eps.eps test1_ADJ.eps 
!gv test1_ADJ.eps

