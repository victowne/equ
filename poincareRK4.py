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
import matplotlib.tri as tri
import h5py
from scipy.spatial import Delaunay

def f(x, y, z, fbr):
    dxdt = fbr(x, y)
    return float(dxdt)

def g(x, y, z, fbz):
    dydt = fbz(x, y)
    return float(dydt)

def e(x, y, z):
    dzdt = 0
    return float(dzdt)

def rk4o(x, y, z, h, fbr, fbz):
    k1x = h*f(x, y, z, fbr)
    k1y = h*g(x, y, z, fbz)
    k1z = h*e(x, y, z)

    k2x = h*f(x + k1x/2.0, y + k1y/2.0, z + k1z/2.0, fbr)
    k2y = h*g(x + k1x/2.0, y + k1y/2.0, z + k1z/2.0, fbz)
    k2z = h*e(x + k1x/2.0, y + k1y/2.0, z + k1z/2.0)

    k3x = h*f(x + k2x/2.0, y + k2y/2.0, z + k2z/2.0, fbr)
    k3y = h*g(x + k2x/2.0, y + k2y/2.0, z + k2z/2.0, fbz)
    k3z = h*e(x + k2x/2.0, y + k2y/2.0, z + k2z/2.0)

    k4x = h*f(x + k3x, y + k3y, z + k3z, fbr)
    k4y = h*g(x + k3x, y + k3y, z + k3z, fbz)
    k4z = h*e(x + k3x, y + k3y, z + k3z)

    x = x + k1x/6.0 + k2x/3.0 + k3x/3.0 + k4x/6.0
    y = y + k1y/6.0 + k2y/3.0 + k3y/3.0 + k4y/6.0
    z = z + k1z/6.0 + k2z/3.0 + k3z/3.0 + k4z/6.0

    return [float(x),float(y),float(z)]

class Mesh():
    def __init__(self, Raxis, Zaxis, Redge, nlayer, fbr, fbz):
        h = .01
        totyList = []
        totxList = []
        totzList = []
        totarclen = []
        numpts = []
        arclenlayer = []
        print('Ready set go')
        for i in np.linspace(Raxis, Redge, nlayer):
            print(i)
            #start point
            yList = []
            xList = []
            zList = []
            arclen = []
            x = i
            y = Zaxis
            z = 0
            xList.append(x)
            yList.append(y)
            zList.append(z)
            t = 1
            changeInTime = h
            #first step
            [x,y,z] = rk4o(xList[t-1], yList[t-1], zList[t-1], h, fbr, fbz)
            xList.append(x)
            yList.append(y)
            zList.append(z)
            temparclen = ((x-xList[t-1])**2 + (y-yList[t-1])**2)**.5
            arclen.append(temparclen)
            t = t + 1
            changeInTime = changeInTime + h

            while changeInTime < 666:

                [x,y,z] = rk4o(xList[t-1], yList[t-1], zList[t-1], h, fbr, fbz)
                if x < xList[t-1] and xList[t-2] < xList[t-1]:
                    break
                xList.append(x)
                yList.append(y)
                zList.append(z)
                temparclen = temparclen + ((x-xList[t-1])**2 + (y-yList[t-1])**2)**.5
                arclen.append(temparclen)
                t = t + 1
                changeInTime = changeInTime + h

            totxList.append(xList[0:-1])
            totyList.append(yList[0:-1])
            totzList.append(zList[0:-1])
            totarclen.append(arclen[0:-1])
            arclenlayer.append(arclen[-1])
            numpts.append(t-1)
        print('find vertex')
        print('layer:1')
        self.Rvtx = []
        self.Zvtx = []
        self.Rvtx.append(Raxis)
        self.Zvtx.append(Zaxis)
        for i in range(1,nlayer):
            print('layer:',i+1)
            dl = arclenlayer[i]/(4*i)
            for j in range(0,4*i):
                temp = list(abs(totarclen[i]-j*dl))
                indexmin = temp.index(min(temp))
                self.Rvtx.append(totxList[i][indexmin])
                self.Zvtx.append(totyList[i][indexmin])

    def pltvtx(self):
        plt.scatter(self.Rvtx,self.Zvtx)
        plt.show()
    
    def crt(self):
        print('create mesh')
        Rvtx = np.array(self.Rvtx)
        Zvtx = np.array(self.Zvtx)
        points = np.vstack((Rvtx,Zvtx)).T
        tri = Delaunay(points)
        meshh5 = h5py.File('mesh.hdf5','w')
        meshh5['points'] = tri.points
        meshh5['triangles'] = tri.simplices
        meshh5['neighbors'] = tri.neighbors
        meshh5.close()
        plt.triplot(Rvtx,Zvtx,tri.simplices.copy())
        plt.plot(Rvtx,Zvtx,'o')
        plt.axis('equal')
        plt.show()

    def itp(self, Psi_bound, Psi_axis, funcpsi, qin, pin, fin, ff, fp):
        print('start interpolation')
        psivtx = []
        qvtx = []
        pvtx = []
        fvtx = []
        Btvtx = []
        Jtvtx = []
        for i in range(0,len(self.Rvtx)):
            temp = float(funcpsi(self.Rvtx[i],self.Zvtx[i]))
            if temp > Psi_bound:
                psivtx.append(Psi_bound)
            elif temp < Psi_axis:
                psivtx.append(Psi_axis)
            else:
                psivtx.append(temp)
        for i in range(0,len(self.Rvtx)):
            qvtx.append(float(qin(psivtx[i])))
            pvtx.append(float(pin(psivtx[i])))
            fvtx.append(float(fin(psivtx[i])))
            Btvtx.append(float(fin(psivtx[i])/self.Rvtx[i]))
            Jtvtx.append(float(self.Rvtx[i]*fp(psivtx[i]) + ff(psivtx[i])/(self.Rvtx[i]*4e-7*np.pi)))
        print('output')
        vtxh5 = h5py.File('equonvtx.hdf5','w')
        vtxh5['R'] = self.Rvtx
        vtxh5['Z'] = self.Zvtx
        vtxh5['psi'] = psivtx
        vtxh5['q'] = qvtx
        vtxh5['p'] = pvtx
        vtxh5['f'] = fvtx
        vtxh5['Bt'] = Btvtx
        vtxh5['Jt'] = Jtvtx
        vtxh5.close()
