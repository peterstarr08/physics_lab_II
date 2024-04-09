import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import tsem
from scipy import odr



_m, _VC, _SG, _temperature = 1, 0.1, 1/50, 0.1

m = 865

#soecific heat J/g/degC
s = 0.33*4.18

diameterData = np.array([111+4*_VC, 111+9*_VC, 111+4*_VC, 111+7*_VC])/1000
widthData = np.array([3+20*_SG, 3+20*_SG, 3+21*_SG, 3+21*_SG])/1000

area = np.pi*np.power(diameterData/2, 2)

a, _a = np.mean(area), tsem(area)
l, _l = np.mean(widthData), tsem(widthData)

print(area, a , _a)


T1, T2 = 83.5, 94.5

coolingData = np.array([
    93.5,
    93,
    92,
    90.5,
    89.5,
    88,
    86,
    87,
    85,
    84.5,
    83.5,
    82.5,
    82,
    81.5,
    81.5,
    81,
    80.5,
    80,
    79.5,
    79,
    78.5,
    78,
    78.5,
    78,
    77.5,
    77.5,
    77,
    76.5,
    76
])
#taken 10 sec intervals
timeData = np.linspace(0, len(coolingData)*10-10, len(coolingData))


ODRdata = odr.RealData(timeData, coolingData, sy=_temperature)
ODRinstantiate = odr.ODR(ODRdata, odr.unilinear)
ODRrun = ODRinstantiate.run()

slope, intercept = ODRrun.beta
_slope, _intercept = ODRrun.sd_beta

xFit = np.linspace(0, 300, 10000)
yFit = xFit*slope + intercept

k = m*s*slope*l/(a*(T2-T1))
print(k)

plt.scatter(timeData, coolingData, s=5)
plt.errorbar(timeData, coolingData, yerr=_temperature, linestyle='')
plt.plot(xFit, yFit)
# plt.show()