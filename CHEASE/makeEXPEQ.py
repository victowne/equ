import numpy as np
import matplotlib.pyplot as plt
nbound = 301 #number of (R,Z)
NWBPS = 2 # =1,only plasma; =2 plasma&wall
NDATA = 2 # =2,uniform resistivity
NPPF1 = 201 # grid number of the profiles
NSTTP = 2 # =1,TT';=2,toroidal current density;=3,parallel current density

#Specify Plasma Boundary#
#Dee Formula
rmax = 0.43
rzero = 1.9
eshape = 1.0 #elongation
xshape = 0.0 #triangularity
tdum = np.linspace(-np.pi,np.pi,nbound);
R = rzero + rmax * (np.cos(tdum) - xshape * np.sin(tdum)**2)
Z = eshape * rmax * np.sin(tdum)
rmax = rmax + 0.1
R2 = rzero + rmax * (np.cos(tdum) - xshape * np.sin(tdum)**2)
Z2 = eshape * rmax * np.sin(tdum)
epsilon = rmax/rzero

#Profiles
s = np.linspace(0,1,NPPF1) # s = sqrt((psi-psi_axis)/(psi_edge-psi_axis))
dpdr = -12*s*(1-s**2)**5
q = 1.33 * (1 + (s/0.595)**8)**.25
ds2 = (s[1] - s[0])*2
dqdr = np.ones(NPPF1)
for i in range(1,NPPF1-1):
    dqdr[i] = (q[i+1] - q[i-1])/ds2
dqdr[0] = 1.33333333*dqdr[1] - 0.33333333*dqdr[2] 
dqdr[NPPF1-1] = 1.33333333*dqdr[NPPF1-2] - 0.33333333*dqdr[NPPF1-3]
dqdr = 1.33 * .25 * (1 + (s/0.595)**8)**-.75 * 8 * (s/0.595)**7/0.595
cur = (2-s/q*dqdr)/q

#Write file
f = open('EXPEQ.test','w')
f.write('  %.16e\n' %epsilon)
f.write('  %.16e\n' %0.)
f.write('  %.16e\n' %0.)
f.write(' %d    %d    %d\n' %(nbound,NWBPS,NDATA))
for i in range(0,nbound):
    f.write('  %.16e %.16e\n' %(R[i],Z[i]))
for i in range(0,nbound):
    f.write('  %.16e %.16e\n' %(R2[i],Z2[i]))
f.write(' %d\n' %NPPF1)
f.write('    %d\n' %NSTTP)
for i in range(0,NPPF1):
    f.write('  %.16e\n' %s[i])
for i in range(0,NPPF1):
    f.write('  %.16e\n' %dpdr[i])
for i in range(0,NPPF1):
    f.write('  %.16e\n' %cur[i])
f.close()
