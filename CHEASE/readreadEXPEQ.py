import numpy as np 
import matplotlib.pyplot as plt
EXP = open('EXPEQ.OUT')
data = EXP.read()
test = data.split()
i=0
a_R0 = float(test[i])
i += 1
Z0_R0 = float(test[i])
i += 1
Pedge = float(test[i])
i += 1
NBPS = int(test[i])
i += 1
NWBPS = int(test[i])
i += 1
NDATA = int(test[i])
i += 1
pbdr = np.zeros((NBPS,NDATA))
wbdr = np.zeros((NBPS,NDATA))
n = 0
for fk in range(i,i+NBPS*2,2):
	pbdr[n,0] = float(test[fk])
	pbdr[n,1] = float(test[fk+1])
	n += 1
	i += 2
n=0
if NWBPS == 2:
	for fk in range(i,i+NBPS*2,2):
		wbdr[n,0] = float(test[fk])
		wbdr[n,1] = float(test[fk+1])
		n +=1
		i += 2
NPPF1 = int(test[i])
i += 1
NSTTP = int(test[i])
i += 1
rgrid = np.zeros(NPPF1)
dpdpsi = np.zeros(NPPF1)
Jtbar = np.zeros(NPPF1)
n=0
for fk in range(i,i+NPPF1):
	rgrid[n] = float(test[fk])
	n +=1
	i += 1
n=0
for fk in range(i,i+NPPF1):
	dpdpsi[n] = float(test[fk])
	n +=1
	i += 1
n=0
for fk in range(i,i+NPPF1):
	Jtbar[n] = float(test[fk])
	n +=1
	i += 1
plt.plot(rgrid,dpdpsi)
plt.figure()
plt.plot(rgrid,Jtbar)
plt.figure()
plt.scatter(pbdr[:,0],pbdr[:,1],s=.5)
plt.scatter(wbdr[:,0],wbdr[:,1],s=.5)
plt.axis('equal')
plt.show()
