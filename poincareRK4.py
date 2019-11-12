import numpy as np
#import matplotlib.pyplot as plt
from reqdsk import funcbr,funcbz
from reqdsk import Raxis,Zaxis,Rmin,Rboxlen,bdr

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

    return [x,y,z]

def main(nlayer):
    totyList = []
    totxList = []
    totzList = []
    totarclen = []
    numpts = []
    arclenlayer = []
    for i in np.linspace(Raxis,bdr[0]-0.01,nlayer):
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
        arclen.append(np.sqrt((x-xList[t-1])**2 + (y-yList[t-1])**2))
        t += 1
        changeInTime += h

        while changeInTime < 60:

            [x,y,z] = rk4o(xList[t-1], yList[t-1], zList[t-1])
            if x < xList[t-1] and xList[t-2] < xList[t-1]:
                break
            xList.append(x)
            yList.append(y)
            zList.append(z)
            arclen.append(np.sqrt((x-xList[t-1])**2 + (y-yList[t-1])**2))
            t += 1
            changeInTime += h

        totxList.append(xList[0:-1])
        totyList.append(yList[0:-1])
        totzList.append(zList[0:-1])
        totarclen.append(arclen[0:-1])
        arclenlayer.append(sum(arclen[0:-1]))
        numpts.append(t-1)
#        plt.scatter(xList[t-3],yList[t-3],marker='x',s=10.)
#        plt.scatter(xList[t-2],yList[t-2],marker='o',s=10.)
#        plt.scatter(xList[t-1],yList[t-1],marker='.',s=10.)
#        plt.show()
    Rvtx = []
    Zvtx = []
    Rvtx.append(Raxis)
    Zvtx.append(Zaxis)
    for i in range(1,nlayer):
        dl = arclenlayer[i]/(4*i)
        for j in range(0,4*i):
           min(totarclen[i]-j*dl) de xia biao 
