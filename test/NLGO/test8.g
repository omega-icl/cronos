set xlabel 't'
set xrange [0:]
#set format x "{/:Italic u}_{%.0f}"; set xtics 1,1,15
set ylabel 'P_{ad}(t)'
set key reverse left Left spacing 2
set bars small

plot 'test8_STA.dat' u 1:(($2+$3)*0.5):(0.5*($3-$2)) tit 'MILP reachability' w yerrorbars lc rgb "blue" lw 1 ps 0
#plot 'test8_STA.dat' u 1:2:3 tit 'MILP reachability' w filledcurves lc rgb "blue"

pause -1 "<ENTER> TO CONTINUE"

set term post eps enh solid color 21
set out 'test8_STA.eps'
rep
!gv test8_STA.eps &

