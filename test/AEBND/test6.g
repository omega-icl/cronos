h=0.05

set xlabel 't'
set xrange [0:]
#set format x "{/:Italic u}_{%.0f}"; set xtics 1,1,15
set ylabel 'u(t)'
set key reverse left Left spacing 2
set bars small

plot 'test6.out' u (($1-1)*h):(($6+$7)*0.5):(0.5*($7-$6)) tit 'diag. block decomp.' w yerrorbars lc rgb "blue" lw 5 ps 0, \
     '' u (($1-1)*h):(($10+$11)*0.5):(0.5*($11-$10)) tit 'recur. block decomp.' w yerrorbars lc rgb "red" lw 5 ps 0
#plot 'linODE_expl.res' u 1:(($4+$5)*0.5):0.2:($5-$4) w boxxyerr fs 2 lw 0, \
#     ''  1:(($2+$3)*0.5):0.2:($3-$2) w boxxyerr fc 1 lw 0
#pause -1
set term post eps enh solid color 21
set out 'test6.eps'
rep
!gv test6.eps &

#set xlabel ''
#set xrange [0.5:15.5]
#set format x "{/:Italic u}_{%.0f}"; set xtics 1,1,15
#set xrange [0.5:15.5]
#set key reverse left Left spacing 2
#set bars small

#plot 'linODE_expl.res' u 1:(($4+$5)*0.5):(0.5*($5-$4)) tit 'recursive block decomp.' w yerrorbars lc 1 lw 15 ps 0, \
#     '' u 1:(($2+$3)*0.5):(0.5*($3-$2)) tit 'no block decomp.' w yerrorbars lc 4 lw 15 ps 0
##plot 'linODE_expl.res' u 1:(($4+$5)*0.5):0.2:($5-$4) w boxxyerr fs 2 lw 0, \
##     ''  1:(($2+$3)*0.5):0.2:($3-$2) w boxxyerr fc 1 lw 0
##pause -1
#set term post eps enh solid color 18
#set out 'linODE_expl.eps'
#rep
#!gv linODE_expl.eps &
