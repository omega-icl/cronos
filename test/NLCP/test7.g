file1 = 'test7.out'
file2 = 'test7b.out'

#set size ratio 1
unset key

#monitorSize=system("xrandr | awk '/\*/{sub(/x/,\",\");print $1; exit}'")
#set macros
#set terminal pngcairo size 1280,720#size @monitorSize
set term post eps enh color 10

################################################################################
#set out 'test7_methanol123.png'
set out 'test7_methanol123.eps'
set multiplot layout 2,3 columnsfirst

set xlabel 'D'
set ylabel 'y_{methanol,1}'
plot file1 u (($1+$2)/2.):(($3+$4)/2.):1:2:3:4 w boxxy lt 1 lc 1, \
       '' u (($1+$2)/2.):(($3+$4)/2.) w p ps .5

set ylabel 'y_{methanol,1}'
plot file2 u 1:2:3 w filledcurves lt 1, \
     '' u 1:2 w l lt 1, '' u 1:3 w l lt 1       

set ylabel 'y_{methanol,2}'
plot file1 u (($1+$2)/2.):(($5+$6)/2.):1:2:5:6 w boxxy lt 1 lc 1, \
       '' u (($1+$2)/2.):(($5+$6)/2.) w p ps .5

set ylabel 'y_{methanol,2}'
plot file2 u 1:4:5 w filledcurves lt 1, \
     '' u 1:4 w l lt 1, '' u 1:5 w l lt 1           

set ylabel 'y_{methanol,3}'
plot file1 u (($1+$2)/2.):(($7+$8)/2.):1:2:7:8 w boxxy lt 1 lc 1, \
       '' u (($1+$2)/2.):(($7+$8)/2.) w p ps .5
       
set ylabel 'y_{methanol,3}'
plot file2 u 1:6:7 w filledcurves lt 1, \
     '' u 1:6 w l lt 1, '' u 1:7 w l lt 1         

unset multiplot

pause -1 "PAUSED"
   
################################################################################
#set out 'test7_methyl-butyrate123.png'
set out 'test7_methyl-butyrate123.eps'
set multiplot layout 2,3 columnsfirst

set xlabel 'D'
set ylabel 'y_{methyl-butyrate,1}'
plot file1 u (($1+$2)/2.):(($9+$10)/2.):1:2:9:10 w boxxy lt 1 lc 1, \
       '' u (($1+$2)/2.):(($9+$10)/2.) w p ps .5

set ylabel 'y_{methyl-butyrate,1}'
plot file2 u 1:8:9 w filledcurves lt 1, \
     '' u 1:8 w l lt 1, '' u 1:9 w l lt 1       

set ylabel 'y_{methyl-butyrate,2}'
plot file1 u (($1+$2)/2.):(($11+$12)/2.):1:2:11:12 w boxxy lt 1 lc 1, \
       '' u (($1+$2)/2.):(($11+$12)/2.) w p ps .5

set ylabel 'y_{methyl-butyrate,2}'
plot file2 u 1:10:11 w filledcurves lt 1, \
     '' u 1:10 w l lt 1, '' u 1:11 w l lt 1    

set ylabel 'y_{methyl-butyrate,3}'
plot file1 u (($1+$2)/2.):(($13+$14)/2.):1:2:13:14 w boxxy lt 1 lc 1, \
       '' u (($1+$2)/2.):(($13+$14)/2.) w p ps .5

set ylabel 'y_{methyl-butyrate,3}'
plot file2 u 1:12:13 w filledcurves lt 1, \
     '' u 1:12 w l lt 1, '' u 1:13 w l lt 1    
          
unset multiplot

pause -1 "PAUSED"

################################################################################
#set out 'test7_toluene123.png'
set out 'test7_toluene123.eps'
set multiplot layout 2,3 columnsfirst

set xlabel 'D'
set ylabel 'y_{toluene,1}'
plot file1 u (($1+$2)/2.):(($15+$16)/2.):1:2:15:16 w boxxy lt 1 lc 1, \
       '' u (($1+$2)/2.):(($15+$16)/2.) w p ps .5

set ylabel 'y_{toluene,1}'
plot file2 u 1:14:15 w filledcurves lt 1, \
     '' u 1:14 w l lt 1, '' u 1:15 w l lt 1       

set ylabel 'y_{toluene,2}'
plot file1 u (($1+$2)/2.):(($17+$18)/2.):1:2:17:18 w boxxy lt 1 lc 1, \
       '' u (($1+$2)/2.):(($17+$18)/2.) w p ps .5

set ylabel 'y_{toluene,2}'
plot file2 u 1:16:17 w filledcurves lt 1, \
     '' u 1:16 w l lt 1, '' u 1:17 w l lt 1 

set ylabel 'y_{toluene,3}'
plot file1 u (($1+$2)/2.):(($19+$20)/2.):1:2:19:20 w boxxy lt 1 lc 1, \
       '' u (($1+$2)/2.):(($19+$20)/2.) w p ps .5

set ylabel 'y_{toluene,3}'
plot file2 u 1:18:19 w filledcurves lt 1, \
     '' u 1:18 w l lt 1, '' u 1:19 w l lt 1 
          
unset multiplot 

pause -1 "PAUSED"

################################################################################

