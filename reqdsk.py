import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy import interpolate
from numpy import gradient as grad
from scipy.interpolate import interp1d
from scipy.integrate import cumtrapz

def parr(f,dR,dZ):
    return grad(f,dR,dZ,edge_order = 2)[1]

def parz(f,dR,dZ):
    return grad(f,dR,dZ,edge_order = 2)[0]  

def not_empty(s):
    return s and s.strip()

with open('EQDSK.OUT','r') as f:
    lines = f.readlines()
#get nr,nz
nr = int(lines[0].split()[-2])
nz = int(lines[0].split()[-1])
#read common arguments
x = []
for line in lines[1:5]:
    x.append(line[0:16])
    x.append(line[16:32])
    x.append(line[32:48])
    x.append(line[48:64])
    x.append(line[64:80])
x = list(filter(not_empty,x))
x = np.array(x,dtype=float)
Rboxlen = x[0]
Zboxlen = x[1]
R0 = x[2]
Rmin = x[3]
Z0 = x[4]
Raxis = x[5]
Zaxis = x[6]
Psi_axis = x[7]
Psi_bound =x[8]
B0 = x[9]
current = x[10]
corp = np.linspace(Psi_axis,Psi_bound,nr)
R = np.ones(nr)
Z = np.ones(nz)
for i in range(0,nr):
    R[i] = Rmin + Rboxlen*i/(nr-1)
    Z[i] = Z0 - Zboxlen*0.5+Zboxlen*i/(nz-1)
#read f
x = []
n = 5
nl = int(nr/5. + 0.8)
for line in lines[n:n+nl]:
    x.append(line[0:16])
    x.append(line[16:32])
    x.append(line[32:48])
    x.append(line[48:64])
    x.append(line[64:80])
x = list(filter(not_empty,x))
f = np.array(x,dtype=float)
#read pressure
n=n+nl
x=[]
for line in lines[n:n+nl]:
    x.append(line[0:16])
    x.append(line[16:32])
    x.append(line[32:48])
    x.append(line[48:64])
    x.append(line[64:80])
x = list(filter(not_empty,x))
p = np.array(x,dtype=float)
#read ffprime
n=n+nl
x=[]
for line in lines[n:n+nl]:
    x.append(line[0:16])
    x.append(line[16:32])
    x.append(line[32:48])
    x.append(line[48:64])
    x.append(line[64:80])
x = list(filter(not_empty,x))
ffprime = np.array(x,dtype=float)
#read pprime
n=n+nl
x=[]
for line in lines[n:n+nl]:
    x.append(line[0:16])
    x.append(line[16:32])
    x.append(line[32:48])
    x.append(line[48:64])
    x.append(line[64:80])
x = list(filter(not_empty,x))
pprime = np.array(x,dtype=float)
#read psi
npsil = int((nr*nz)/5. + 0.8)
n=n+nl
x=[]
for line in lines[n:n+npsil]:
    x.append(line[0:16])
    x.append(line[16:32])
    x.append(line[32:48])
    x.append(line[48:64])
    x.append(line[64:80])
x = list(filter(not_empty,x))
psi = np.zeros((nr,nz))
ntest = 0
for i in range(0,nz):
    for j in range(0,nr):
        psi[i,j]=float(x[ntest])
        ntest  = ntest + 1
#read q
n = n + npsil
x = []
for line in lines[n:n+nl]:
    x.append(line[0:16])
    x.append(line[16:32])
    x.append(line[32:48])
    x.append(line[48:64])
    x.append(line[64:80])
x = list(filter(not_empty,x))
q = np.array(x,dtype=float)
#read plasma boundary
n = n + nl
nbdr = int(lines[n].split()[0])
nlmt = int(lines[n].split()[1])
n = n + 1
nl2 = int((nbdr + nlmt)*.4 + 0.8)
x = []
for line in lines[n:n+nl2]:
    x.append(line[0:16])
    x.append(line[16:32])
    x.append(line[32:48])
    x.append(line[48:64])
    x.append(line[64:80])
x = list(filter(not_empty,x))
pts = np.array(x,dtype=float).reshape(-1,2)
bdr = pts[0:nbdr,0]
bdz = pts[0:nbdr,1]
lmr = pts[nbdr:nbdr+nlmt,0]
lmz = pts[nbdr:nbdr+nlmt,1]
#start_interpolation
dR = R[1]-R[0]
dZ = Z[1]-Z[0]
Br =  parz(psi,dZ,dR)/R
Bz = -parr(psi,dZ,dR)/R
funcpsi = interpolate.RectBivariateSpline(R,Z,psi.T)
funcbr = interpolate.RectBivariateSpline(R,Z,Br.T)
funcbz = interpolate.RectBivariateSpline(R,Z,Bz.T)
qin = interp1d(corp,q,kind = 'cubic')
pin = interp1d(corp,p,kind = 'cubic')
fin = interp1d(corp,f,kind = 'cubic')
ff = interp1d(corp,ffprime,kind = 'cubic')
fp = interp1d(corp,pprime,kind = 'cubic')
psi_p = np.linspace(0,1,nr)
psit = cumtrapz(q,psi_p,initial=0)
psi_norm = psit - min(psit)
psi_norm = psi_norm/psi_norm[nr-1]
r = np.sqrt(psi_norm)
#eqdsk -> hdf5
eqh5 = h5py.File('equlibrium.hdf5','w')
eqh5.attrs['B0'] = B0
eqh5.attrs['current'] = current
eqh5.attrs['Raxis'] = Raxis
eqh5.attrs['Zaxis'] = Zaxis
eqh5['r'] = r
eqh5['psi_1d'] = corp
eqh5['psi_norm'] = psi_norm 
eqh5['f'] = f
eqh5['p'] = p
eqh5['q'] = q
eqh5['ffprime'] = ffprime
eqh5['pprime'] = pprime
eqh5['R'] = R
eqh5['Z'] = Z
eqh5['psi'] = psi
eqh5['Br'] = Br
eqh5['Bz'] = Bz
eqh5['boundary'] = [bdr,bdz]
eqh5['limiter'] = [lmr,lmz]
eqh5.close()
