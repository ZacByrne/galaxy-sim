from pylab import *
import numpy as np
data = loadtxt("dwarf")
r = data[:,0]
r = r /(6.171*10.0**19.0)
vc = data[:,3]
sigma2 = data[:,4]
plot(np.log10((r)), np.log10(vc), label="Numerical")
plot(np.log10(r),np.log10(sigma2), label ="Analytic")
#plot(np.log10(r), (np.log10(sigma2)-np.log10(vc)), label = "Difference between black hole and no black hole")
xlabel('(Log10(R/R_c))       ( R_c = 2 Kpc)')
ylabel('log10 Velocity dispersion (log10(m/s))')
title('Plummer model velocity dispersion')
xlim([-2,1.5])
legend(loc = 3)
show()

plot(np.log10(r), (100*(sigma2 -vc))/(sigma2), label = "%(analytic - numeric)/analytic")
xlabel(' (Log10(R/R_c)) R_c = 5 Kpc')
ylabel('Percentage difference')
title('Plummer model percentage difference between analytic and numeric')
xlim([-2,1.5])
ylim([-0.01,0.01])
legend(loc = 2)
show()
