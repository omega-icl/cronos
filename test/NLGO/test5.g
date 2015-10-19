tau = 5.5e-3
kr  = 1.4e-4
kd  = 5.0e-6
K    = kd/kr #3.6e-2

alpha = 0.0160
beta  = 0.492
kappa = 0.469

thetaLL = 0.082
thetaHL = 0.018

sigma(theta) = exp(kappa*log(theta))
mu(I,theta)  = alpha*I/(1.+tau*beta*sigma(theta)*I + K*tau*(beta*sigma(theta)*I)**2)
muopt(I) = alpha*sqrt(K*tau)/(tau+2*sqrt(K*tau))*I

set xrange [0:1400]
set xlabel 'I, muE/m2/s'
set yrange [0:7]
set ylabel 'mu, gC/gchl/h'
set grid

plot mu(x,thetaLL) tit 'LL' w l lt 1 lc 1 lw 3, mu(x,thetaHL) tit 'HL' w l lt 1 lc 3 lw 3, muopt(x) tit '' w l lt 3 lc -1 lw 3

pause -1
