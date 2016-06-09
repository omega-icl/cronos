file1 = 'mp-NLP1.out'
file2 = 'mp-NLP2.out'
file3 = 'mp-NLP3.out'
file4 = 'mp-NLP4.out'

set size ratio 1
set key out reverse Left spacing 2

set xlabel 'x1'
set ylabel 'x2'
plot file1 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 tit 'no act.' w boxxy lt 1 lc 1 fs transparent solid 0.3 nobo, \
     file2 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 tit 'g1, g2 act.' w boxxy lt 1 lc 2 fs transparent solid 0.3 nobo, \
     file3 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 tit 'g1 act.' w boxxy lt 1 lc 3 fs transparent solid 0.3 nobo, \
     file4 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 tit 'g2 act.' w boxxy lt 1 lc 4 fs transparent solid 0.3 nobo

#pause -1
set term pdfcairo enh color font "Arial,12"
set out 'mp-NLP.pdf'
rep
set term x11
!pdfcrop mp-NLP.pdf
!mv mp-NLP-crop.pdf mp-NLP.pdf
!gv mp-NLP.pdf
#set term post eps enh color 18
#set out 'mp-NLP.eps'
#rep
#set term x11
#!ps2eps -B -f -l mp-NLP.eps
#!mv mp-NLP.eps.eps mp-NLP.eps


