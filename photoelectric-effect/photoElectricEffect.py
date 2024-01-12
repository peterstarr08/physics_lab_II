import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

e = 1.6e-19
c = 299792458.

#note x is frequency
def stoppingPotential(x, h, phi):
    return (h/e)*x - phi

wavelengthData = np.array([635, 570, 540, 500, 460])
stoppingPotentialData, potentialUncertainty = np.array([-0.42, -0.63, -0.82, -0.94, -1.20]), 0.01

params, covar = curve_fit(stoppingPotential, c*(10**9)/wavelengthData, np.absolute(stoppingPotentialData), sigma=np.full(len(stoppingPotentialData), potentialUncertainty))
standardDeviations = np.sqrt(np.diag(covar))

print("h: "+str(params[0]) + " +/- "+ str(standardDeviations[0]))
print("work function: "+str(params[1]) + " +/- "+ str(standardDeviations[1]) + " eV")

plt.scatter(c*(10**9)/wavelengthData, np.absolute(stoppingPotentialData), s=5)
plt.errorbar(c*(10**9)/wavelengthData, np.absolute(stoppingPotentialData), potentialUncertainty, linestyle='')

xFit = np.linspace(4.70e14,6.55e14, 10000)
yFit = (params[0]/e)*xFit - params[1]

plt.plot(xFit, yFit, label='Regression curve')
# plt.legend(loc='top left')

plt.show()
