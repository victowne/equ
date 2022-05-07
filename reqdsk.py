'''
    Copyright (C) <2020> <Author: Weikang Tang>
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
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

class gfile:

    def __init__(self,filename):
        try:
            with open(filename,'r') as f:
                lines = f.readlines()
            print('open gfile successfully')
        except:
            print('filename not exist')
            return
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
        self.Rboxlen = x[0]
        self.Zboxlen = x[1]
        self.R0 = x[2]
        self.Rmin = x[3]
        self.Z0 = x[4]
        self.Raxis = x[5]
        self.Zaxis = x[6]
        self.Psi_axis = x[7]
        self.Psi_bound =x[8]
        self.B0 = x[9]
        self.current = x[10]
        corp = np.linspace(self.Psi_axis,self.Psi_bound,nr)
        self.psip1d = corp
        R = np.ones(nr)
        Z = np.ones(nz)
        for i in range(0,nr):
            R[i] = self.Rmin + self.Rboxlen*i/(nr-1)
            Z[i] = self.Z0 - self.Zboxlen*0.5+self.Zboxlen*i/(nz-1)
        self.R = R
        self.Z = Z
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
        self.f = np.array(x,dtype=float)
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
        self.p = np.array(x,dtype=float)
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
        self.ffprime = np.array(x,dtype=float)
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
        self.pprime = np.array(x,dtype=float)
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
        self.psi = np.zeros((nr,nz))
        ntest = 0
        for i in range(0,nz):
            for j in range(0,nr):
                self.psi[i,j]=float(x[ntest])
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
        self.q = np.array(x,dtype=float)
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
        self.bdr = pts[0:nbdr,0]
        self.bdz = pts[0:nbdr,1]
        self.lmr = pts[nbdr:nbdr+nlmt,0]
        self.lmz = pts[nbdr:nbdr+nlmt,1]
        print('read gfile finished')
        #start_interpolation
        dR = R[1]-R[0]
        dZ = Z[1]-Z[0]
        self.Br =  parz(self.psi,dZ,dR)/R
        self.Bz = -parr(self.psi,dZ,dR)/R
        self.funcpsi = interpolate.RectBivariateSpline(R,Z,self.psi.T)
        self.funcbr = interpolate.RectBivariateSpline(R,Z,self.Br.T)
        self.funcbz = interpolate.RectBivariateSpline(R,Z,self.Bz.T)
        self.qin = interp1d(corp,self.q,kind = 'cubic')
        self.pin = interp1d(corp,self.p,kind = 'cubic')
        self.fin = interp1d(corp,self.f,kind = 'cubic')
        self.ff = interp1d(corp,self.ffprime,kind = 'cubic')
        self.fp = interp1d(corp,self.pprime,kind = 'cubic')
        psi_p = np.linspace(0,1,nr)
        psit = cumtrapz(self.q,psi_p,initial=0)
        psi_norm = psit - min(psit)
        psi_norm = psi_norm/psi_norm[nr-1]
        self.r = np.sqrt(psi_norm)
        print('interpolation finished')

    def g2h5(self):
        #eqdsk -> hdf5
        eqh5 = h5py.File('equlibrium.hdf5','w')
        eqh5.attrs['B0'] = self.B0
        eqh5.attrs['current'] = self.current
        eqh5.attrs['Raxis'] = self.Raxis
        eqh5.attrs['Zaxis'] = self.Zaxis
        eqh5['r'] = self.r
        eqh5['psi_1d'] = self.psip1d
        eqh5['f'] = self.f
        eqh5['p'] = self.p
        eqh5['q'] = self.q
        eqh5['ffprime'] = self.ffprime
        eqh5['pprime'] = self.pprime
        eqh5['R'] = self.R
        eqh5['Z'] = self.Z
        eqh5['psi'] = self.psi
        eqh5['Br'] = self.Br
        eqh5['Bz'] = self.Bz
        eqh5['boundary'] = [self.bdr,self.bdz]
        eqh5['limiter'] = [self.lmr,self.lmz]
        eqh5.close()

if __name__ == '__main__':
    gfile(sys.argv[1]).g2h5
