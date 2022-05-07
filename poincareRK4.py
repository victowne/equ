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
from reqdsk import funcbr,funcbz,funcpsi
from reqdsk import Raxis,Zaxis,Rmin,Rboxlen,bdr,Psi_axis,Psi_bound
from reqdsk import qin,pin,fin,ff,fp

h = .01

def f(x,y,z):
    dxdt = funcbr(x,y)
    return dxdt

def g(x,y,z):
    dydt = funcbz(x,y)
    return dydt

def e(x,y,z):
    dzdt = 0
    return dzdt

def rk4o(x, y, z):
    global h
    k1x = h*f(x, y, z)
    k1y = h*g(x, y, z)
    k1z = h*e(x, y, z)

    k2x = h*f(x + k1x/2.0, y + k1y/2.0, z + k1z/2.0)
    k2y = h*g(x + k1x/2.0, y + k1y/2.0, z + k1z/2.0)
    k2z = h*e(x + k1x/2.0, y + k1y/2.0, z + k1z/2.0)

    k3x = h*f(x + k2x/2.0, y + k2y/2.0, z + k2z/2.0)
    k3y = h*g(x + k2x/2.0, y + k2y/2.0, z + k2z/2.0)
    k3z = h*e(x + k2x/2.0, y + k2y/2.0, z + k2z/2.0)

    k4x = h*f(x + k3x, y + k3y, z + k3z)
    k4y = h*g(x + k3x, y + k3y, z + k3z)
    k4z = h*e(x + k3x, y + k3y, z + k3z)

    x = x + k1x/6.0 + k2x/3.0 + k3x/3.0 + k4x/6.0
    y = y + k1y/6.0 + k2y/3.0 + k3y/3.0 + k4y/6.0
    z = z + k1z/6.0 + k2z/3.0 + k3z/3.0 + k4z/6.0

    return [float(x),float(y),float(z)]

def main(nlayer):
    totyList = []
    totxList = []
    totzList = []
    totarclen = []
    numpts = []
    arclenlayer = []
    print('Ready set go')
    for i in np.linspace(Raxis,bdr[0],nlayer):
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
        [x,y,z] = rk4o(xList[t-1], yList[t-1], zList[t-1])
        xList.append(x)
        yList.append(y)
        zList.append(z)
        temparclen = ((x-xList[t-1])**2 + (y-yList[t-1])**2)**.5
        arclen.append(temparclen)
        t = t + 1
        changeInTime = changeInTime + h

        while changeInTime < 666:

            [x,y,z] = rk4o(xList[t-1], yList[t-1], zList[t-1])
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
        arclenlayer.append(arclen[-2])
        numpts.append(t-1)
    print('find vertex')
    Rvtx = []
    Zvtx = []
    Rvtx.append(Raxis)
    Zvtx.append(Zaxis)
    for i in range(1,nlayer):
        print(i+1)
        dl = arclenlayer[i]/(4*i)
        for j in range(0,4*i):
            temp = list(abs(totarclen[i]-j*dl))
            indexmin = temp.index(min(temp))
            Rvtx.append(totxList[i][indexmin])
            Zvtx.append(totyList[i][indexmin])
    print('start interpolation')
    psivtx = []
    qvtx = []
    pvtx = []
    fvtx = []
    Btvtx = []
    Jtvtx = []
    for i in range(0,len(Rvtx)):
        temp = float(funcpsi(Rvtx[i],Zvtx[i]))
        if temp > Psi_bound:
            psivtx.append(Psi_bound)
        elif temp < Psi_axis:
            psivtx.append(Psi_axis)
        else:
            psivtx.append(temp)
    for i in range(0,len(Rvtx)):
        qvtx.append(float(qin(psivtx[i])))
        pvtx.append(float(pin(psivtx[i])))
        fvtx.append(float(fin(psivtx[i])))
        Btvtx.append(float(fin(psivtx[i])/Rvtx[i]))
        Jtvtx.append(float(Rvtx[i]*fp(psivtx[i]) + ff(psivtx[i])/(Rvtx[i]*4e-7*np.pi)))
    print('output')
    vtxh5 = h5py.File('equonvtx.hdf5','w')
    vtxh5['R'] = Rvtx
    vtxh5['Z'] = Zvtx
    vtxh5['psi'] = psivtx
    vtxh5['q'] = qvtx
    vtxh5['p'] = pvtx
    vtxh5['f'] = fvtx
    vtxh5['Bt'] = Btvtx
    vtxh5['Jt'] = Jtvtx
    vtxh5.close()
    Rvtx = np.array(Rvtx)
    Zvtx = np.array(Zvtx)
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
