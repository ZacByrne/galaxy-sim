from pylab import *
import numpy as np
data = loadtxt("analy")

r = data[:,0]
r = r /(1.543*10.0**20.0)
vc = data[:,1]
sigma2 = data[:,2]

data = loadtxt("ai1")
bh = data[:,1]

data = loadtxt("ai05")
bhha = data[:,1]

data = loadtxt("05ai")
bhh = data[:,1]

data = loadtxt("1ai")
aii = data[:,1]

plot(np.log10((r)), np.log10(vc), label="Numerical beta = 0")
plot(np.log10(r),np.log10(bh), label ="Numerical beta = 1")
plot(np.log10(r),np.log10(bhha), label ="Numerical beta = 0.5")
plot(np.log10(r),np.log10(bhh), label ="Numerical beta = -0.5")
plot(np.log10(r),np.log10(aii), label ="Numerical beta = -1")
#plot(np.log10(r), (np.log10(sigma2)-np.log10(vc)), label = "Difference between black hole and no black hole")
xlabel('log10 Radius (Log10(R_c)) R_c = 5 Kpc')
ylabel('log10 Velocity dispersion (log10(m/s))')
title('Plummer model velocity dispersion with anisotropy')
xlim([-2,2])
legend(loc = 3)
show()

