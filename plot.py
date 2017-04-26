from pylab import *
data = loadtxt("test")
r = data[:,0]
vc = data[:,1]
#sigma2 = data[:,2]
plot(r, vc)
#plot(r,sigma2)
show()
