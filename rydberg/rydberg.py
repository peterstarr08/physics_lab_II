import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import tstd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

LC = 1/60
GRATING_SPACE = 1/15000 * 0.0254
R = 1.097e7

def regression(x, R):
    return R*x

m = np.array([1,1,1,1, 2])
n_i = np.array([6, 5,4 ,3,4])

vernier_1_tL =np.radians( 360- np.array([347.5+2*LC, 346.5+11*LC, 346 +0*LC, 341+7*LC, 331+10*LC]))
vernier_1_tR = np.radians(np.array([12.5+1*LC, 13+21*LC, 14+ 0*LC, 19+6*LC, 29+27*LC]))

vernier_2_tL = np.radians(180-np.array([167.5+15*LC, 166.5+17*LC, 166+11*LC, 161 + 7*LC, 151+26*LC]))
vernier_2_tR = np.radians(np.array([192.5 + 1*LC, 193+22*LC, 194 + 0*LC, 199 + 4*LC, 206.5 + 21*LC])-180)

rwavelength_v1_l = 1/(GRATING_SPACE*np.sin(vernier_1_tL)/m)
rwavelength_v1_r = 1/(GRATING_SPACE*np.sin(vernier_1_tR)/m)
rwavelength_v2_l = 1/(GRATING_SPACE*np.sin(vernier_2_tL)/m)
rwavelength_v2_r = 1/(GRATING_SPACE*np.sin(vernier_2_tR)/m)


_wavelength = np.full(len(n_i), .0)
deviationError = np.full(len(_wavelength), .0)
for i in range(0, len(_wavelength)):
    _wavelength[i] = (rwavelength_v1_l[i] + rwavelength_v1_r[i] + rwavelength_v2_l[i] + rwavelength_v2_r[i])/4
    deviationError[i] = tstd([rwavelength_v1_l[i] , rwavelength_v1_r[i] , rwavelength_v2_l[i] , rwavelength_v2_r[i]])

stateChange = (1/4 - 1/(np.power(n_i, 2)))

params, covar = curve_fit(regression, stateChange, _wavelength, sigma=deviationError)

print("R: ",params[0], "+/-", np.sqrt(np.diag(covar))[0])

print("Percent error: ", abs(R-params[0])/R*100, "%")

plt.title(r"$\frac{1}{\lambda}\text{ vs } \left( \frac{1}{4}-\frac{1}{n_i^{2}} \right) $")
plt.xlabel(r"$\left( \frac{1}{4}-\frac{1}{n_i^{2}} \right) $")
plt.ylabel(r"$\frac{1}{\lambda}(m^{-1})$")
plt.scatter(stateChange, _wavelength, s=5)
plt.errorbar(stateChange, _wavelength, yerr=deviationError, linestyle='')

xFit = np.linspace(0.13, 0.23, 1000)
yFit = params[0]*xFit

plt.plot(xFit, yFit, label='Regression Curve')
plt.legend(loc='upper left')

plt.show()