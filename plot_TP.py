import numpy as np
import matplotlib.pylab as plt
import pandas as pd

data1 = np.loadtxt('FMS_RC_ic.out')
Pic = data1[:,1]
Tic = data1[:,2]

data2 = np.loadtxt('FMS_RC_pp.out')
Ppp = data2[:,1]
Tpp = data2[:,2]
dT_rad = data2[:,3]
dT_conv = data2[:,4]
tauV = data2[:,5]
tauIR = data2[:,6]


fig = plt.figure()

plt.plot(Tic,Pic/1e5,ls='dashed',lw=3,label='Initial Conditions',c='orange')
plt.plot(Tpp,Ppp/1e5,lw=3,label='Numerical Result',c='blue')

plt.ylabel('Pressure [bar]')
plt.xlabel('Temperature [K]')
plt.legend()

plt.yscale('log')
plt.gca().invert_yaxis()

fig = plt.figure()

plt.plot(dT_rad,Ppp/1e5,lw=3,label='dT_rad',c='red')
plt.plot(dT_conv,Ppp/1e5,lw=3,label='dT_conv',c='blue')
plt.ylabel('Pressure [bar]')
plt.xlabel('dT [K s$^{-1}$]')

plt.legend()
plt.yscale('log')
plt.gca().invert_yaxis()

fig = plt.figure()

plt.plot(tauV,Ppp/1e5,lw=3,label='tauV',c='blue')
plt.plot(tauIR,Ppp/1e5,lw=3,label='tauIR',c='red')

plt.ylabel('Pressure [bar]')
plt.xlabel('Optical Depth')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.gca().invert_yaxis()



plt.show()
