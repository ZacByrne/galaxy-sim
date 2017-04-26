from pylab import *
import numpy as np
data = loadtxt("analy")
r = data[:,0]
r = r /(1.543*10.0**20.0)
vc = data[:,1]
sigma2 = data[:,2]
plot(np.log10((r)), np.log10(vc), label="Numerical")
plot(np.log10(r),np.log10(sigma2), label ="Analytic")
#plot(np.log10(r), (np.log10(sigma2)-np.log10(vc)), label = "Difference between black hole and no black hole")
xlabel('log10 Radius (Log10(R_c)) R_c = 5 Kpc')
ylabel('log10 Velocity dispersion (log10(m/s))')
title('Plummer model velocity dispersion')
legend(loc = 3)
show()

plot(np.log10(r), (100*(sigma2 -vc))/(sigma2), label = "%(analytic - numeric)/analytic")
xlabel('log10 Radius (Log10(R_c)) R_c = 5 Kpc')
ylabel('Percentage difference')
title('Plummer model percentage difference between analytic and numeric')
xlim([-2,2])
ylim([-0.01,0.01])
legend(loc = 2)
show()
