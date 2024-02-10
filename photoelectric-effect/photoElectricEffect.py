import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import odr

e = 1.6e-19
c = 299792458.0

#note x is frequency
def stoppingPotential(x, h, phi):
    return (h/e)*x - phi

wavelengthData = np.array([635, 570, 540, 500, 460])
stoppingPotentialData, potentialUncertainty = np.array([-0.42, -0.63, -0.82, -0.94, -1.20]), 0.01

distanceData, distanceUncertainity = np.array(range(40, 16, -2), dtype=float), 0.01
currentData, currenUncertainity = np.array([0.294, 0.316, 0.349, 0.399, 0.407, 0.440, 0.490,0.569, 0.673, 0.814, 1.020, 1.340]), 0.001


errorDistance = []
for val in distanceData:
    errorDistance.append(distanceUncertainity/(val*val*val))
inverseSquareLawData = odr.RealData(1/np.power(distanceData,2), currentData, errorDistance, currenUncertainity)
odrInstantiate = odr.ODR(inverseSquareLawData, odr.unilinear)
val = odrInstantiate.run()

m, intercept = val.beta
dm,dc = val.sd_beta

plt.scatter(1/np.power(distanceData,2), currentData, s=5)
plt.errorbar(1/np.power(distanceData,2), currentData,yerr= currenUncertainity, xerr=errorDistance, linestyle='',  capsize=2)
plt.title(r"$\text{Inverse Square Law}$")
plt.ylabel(r"$\text{Current(}\mu\text{A})$")
plt.xlabel(r"$\frac{1}{r^{2}}\text{ in }cm^{-2}$")
xFit = np.linspace(0.0005, 0.0035, 10000)
yFit = m*xFit + intercept

plt.plot(xFit, yFit, label = "regression curve")
plt.savefig('Inverse_square_law.png', dpi=300)

plt.close()
plt.clf()

params, covar = curve_fit(stoppingPotential, c*(10**9)/wavelengthData, np.absolute(stoppingPotentialData), sigma=np.full(len(stoppingPotentialData), potentialUncertainty))
standardDeviations = np.sqrt(np.diag(covar))

print("h: "+str(params[0]) + " +/- "+ str(standardDeviations[0]))
print("work function: "+str(params[1]) + " +/- "+ str(standardDeviations[1]) + " eV")
plt.scatter(c*(10**9)/wavelengthData, np.absolute(stoppingPotentialData), s=5)
plt.errorbar(c*(10**9)/wavelengthData, np.absolute(stoppingPotentialData), potentialUncertainty, linestyle='', capsize=2)

xFit = np.linspace(4.70e14,6.55e14, 10000)
yFit = (params[0]/e)*xFit - params[1]

plt.plot(xFit, yFit, label='Regression curve')
plt.title(r"$\text{Plancks's constant - Stoping Potenial(V) vs frequency(}\nu)$")
plt.xlabel(r"$\nu(\text{ in } s^{-1})$")
plt.ylabel(r"$\text{Retarding potential (eV) }$")
plt.legend(loc='upper left')

plt.savefig('plancks_constant.png', dpi=300)



