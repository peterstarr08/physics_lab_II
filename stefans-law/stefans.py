import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import odr

_V, _I = 0.01, 0.01

#Voltage Current Temperature
DATA = np.array([
    [0.61, 0.07, 21.8],
    [0.73, 0.08, 22.0],
    [0.84, 0.09, 22.1],
    [1.01, 0.10, 22.3],
    [1.15, 0.11, 22.5],
    [1.24, 0.11, 22.7],
    [1.34, 0.12, 22.9],
    [1.46, 0.12, 23.1],
    [1.53, 0.12, 23.3],
    [1.62, 0.13, 23.5],
    [1.72, 0.14, 23.8],
    [1.94, 0.15, 24.0],
    [2.09, 0.15, 24.2],
    [2.20, 0.16, 24.5],
    [2.35, 0.16, 24.7],
    [2.48, 0.17, 25.0],
    [2.55, 0.17, 25.5],
    [2.64, 0.17, 25.7],
    [2.85, 0.18, 25.9],
    [2.96, 0.19, 26.3],
    [3.11, 0.19, 26.6],
    [3.27, 0.20, 27.0],
    [3.43, 0.20, 27.2],
    [3.77, 0.22, 28.1]
    
])
#Converts to kelvin
DATA[:,2] = DATA[:,2] + 273


alpha = 5.21e-3
beta = 7.2e-7

T_G = 800 #Kelvin

R_G = DATA[0,0]/DATA[0,1]
_R_G = R_G*math.sqrt((_V/DATA[0,0])**2 + (_I/DATA[0,1])**2)

R_0 = (R_G)/(1+alpha*T_G+beta*(T_G**2))#Ohms
_R_0 = abs(_R_G/(1+alpha*T_G+beta*(T_G**2)))

R = DATA[:, 0]/DATA[:,1]
_R  = R*np.sqrt(np.power(_V/DATA[:,0],2)+np.power(_I/DATA[:,1],2))

# print(R, _R)

delta = np.sqrt(alpha**2 - 4*beta*(1-R/R_0))
T = (-1*alpha + delta)/(2*beta)
_T = np.sqrt(np.power((1/4/beta/(delta)*(4*beta/R_0))*_R,2) + np.power((1/4/beta/(delta)*(-4*beta*R/(R_0**2)))*_R_0,2))
# _T =abs((1/4/beta/delta*(4*beta/R_0))*_R) + abs((1/4/beta/delta*(-4*beta*R/(R_0**2)))*_R_0)


P = DATA[:, 0]*DATA[:,1]
_P = P*np.sqrt(np.power(_V/DATA[:,0],2)+np.power(_I/DATA[:,1],2))

xerr = abs(_T/T)
yerr = abs(_P/P)

ODRdata = odr.RealData(np.log(T), np.log(P), sx=xerr, sy=yerr)
ODRinstantiate = odr.ODR(ODRdata, odr.unilinear)
ODRrun = ODRinstantiate.run()

slope, intercept = ODRrun.beta
_slope, _intercept = ODRrun.sd_beta

xFit = np.linspace(np.log(T)[0], np.log(T)[len(T)-1], 1000)
yFit = slope*xFit + intercept

print(slope, '+/-', _slope)
stephan_e = math.exp(intercept)
_stephan_e = math.exp(_intercept)

print(stephan_e/(0.25e-6)/0.950)

plt.errorbar(np.log(T), np.log(P), xerr=xerr, yerr=yerr, linestyle='', capsize=2)
plt.scatter(np.log(T), np.log(P))

plt.plot(xFit, yFit)

# plt.errorbar(T,P, xerr=_T, yerr=_P, linestyle='')
plt.show()