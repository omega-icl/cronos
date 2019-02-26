file1  = 'test4_bnd.out'
file2  = 'test4_inn.out'
NP = 4

set term post eps enh color 8
set out 'test4.eps'

#set size ratio 1
set style fill transparent solid 0.5 noborder
unset key

set multiplot layout NP-1,NP-1 columnsfirst margins 0.1,0.9,0.1,0.9 spacing 0.1
do for [i = 1:NP-1]{
  do for [j = 2:i]{
    set multiplot next
  }
  do for [j = i+1:NP]{
    set xlabel sprintf("p%d",i)
    set ylabel sprintf("p%d",j)
    plot file1 u ((column(2*i-1)+column(2*i))/2.):((column(2*j-1)+column(2*j))/2.):(column(2*i-1)):(column(2*i)):(column(2*j-1)):(column(2*j)) w boxxy lt 1 lc 2, \
         file2 u ((column(2*i-1)+column(2*i))/2.):((column(2*j-1)+column(2*j))/2.):(column(2*i-1)):(column(2*i)):(column(2*j-1)):(column(2*j)) w boxxy lt 1 lc 1
  }
}
unset multiplot

!ps2eps -B -f -l test4.eps
!mv test4.eps.eps test4.eps
!gv test4.eps

