from pylab import *
import numpy as np
data = loadtxt("test")

r = data[:,0]
r = r /(1.543*10.0**20.0)
vc = data[:,3]
sigma2 = data[:,2]

data = loadtxt("bhb")
bh = data[:,3]

data = loadtxt("bha")
bhha = data[:,3]

data = loadtxt("bhten")
bhh = data[:,3]

plot(np.log10((r)), np.log10(vc), label="Numerical without BH")
plot(np.log10(r),np.log10(bh), label ="Numerical with 10% BH")
plot(np.log10(r),np.log10(bhha), label ="Numerical with 1% BH")
plot(np.log10(r),np.log10(bhh), label ="Numerical with .1% BH")
#plot(np.log10(r), (np.log10(sigma2)-np.log10(vc)), label = "Difference between black hole and no black hole")
xlabel('log10 Radius (Log10(R_c)) R_c = 5 Kpc')
ylabel('log10 Velocity dispersion (log10(m/s))')
title('Plummer model velocity dispersion with BH')
xlim([-2,2])
legend(loc = 3)
show()

