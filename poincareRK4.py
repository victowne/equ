import numpy as np
import matplotlib.pyplot as plt

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


for i in np.arange(3.1,3.2,0.1):
    yList = []
    xList = []
    zList = []
    #start point
    x = i
    y = 0
    z = 0
    h = .01
    xList.append(x)
    yList.append(y)
    zList.append(z)
    t = 1
    changeInTime = h

    while changeInTime < 60:

        [x,y,z] = rk4o(xList[t-1], yList[t-1], zList[t-1])
        xList.append(x)
        yList.append(y)
        zList.append(z)
        if 1 < changeInTime:
            if x < xList[t-1] and xList[t-2] < xList[t-1]:
                break

        t += 1
        changeInTime += h
    plt.scatter(xList,yList,s=.1)
plt.axis('equal')
plt.show()
# f = open("location.dat", "w")
# i = 0
# while i < 100:
#       f.write(str(xList[i]) + "," + str(yList[i])+ "," + str(zList[i])+"\n")
#       i = i + 1
# f.close()
