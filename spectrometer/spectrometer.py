import matplotlib.pyplot as plt
import  numpy as np
import math
from tabulate import tabulate
from scipy import odr


LC = 1/60
d = 1/15000 * 0.0254

colorsData = ['Violet', 'Green', 'Yellow','Red']
colorsWavelength = (np.array([407.8, 546.1, 577, 623]))*1e-9

wavelengthData = np.array([
    [347+28*LC, 167.5+6*LC, 12.5+0*LC, 192.5+6*LC],
    [344+12*LC, 164+20*LC, 15.5+15/60, 195.5+18/60],
    [343+10*LC, 163+25*LC, 16.5+11*LC, 196.5+17*LC],
    [342 + 0*LC, 162+6*LC, 18+28*LC, 198+0*LC]
])

prismAngleData = np.array([122+19*LC, 241.5+13*LC])
centralMaxima = 179+10*LC

prismData = np.array([
    139.5+16*LC, 140.5+2*LC, 140.5 + 9*LC, 140.5+17*LC
])


#Normalizing data
wavelengthData[:,0] = np.radians(360-wavelengthData[:, 0])
wavelengthData[:,1] = np.radians(180 - wavelengthData[:,1])
wavelengthData[:,2] = np.radians(wavelengthData[:,2])
wavelengthData[:,3] = np.radians(wavelengthData[:,3]-180)

#Calculation
angleOfPrism = math.radians(((180-prismAngleData[0])+(prismAngleData[1]-180))/2)
deviationAngle = np.radians(centralMaxima - prismData)

uncertaintyAngle = LC

wavelengthMatrix = (d*np.sin(wavelengthData))
inverseSquareWavelengthMatrix = 1/np.power(wavelengthMatrix, 2)

meanWavelengths = np.mean(wavelengthMatrix, axis=1)
meanInverseSquareWavelengths = np.mean(inverseSquareWavelengthMatrix, axis=1)

standardDeviationWavelengths = np.std(wavelengthMatrix, axis=1)
standardDeviationInverseSquareWavelengths = np.std(inverseSquareWavelengthMatrix, axis=1)


refractiveIndex = np.sin((deviationAngle+angleOfPrism)/2)/math.sin(angleOfPrism/2)
uncertaintyRefractiveIndex = uncertaintyAngle*np.sqrt((1/math.sin(angleOfPrism/2)*np.cos(deviationAngle/2 + angleOfPrism/2)/2)**2  +  ((math.sin(angleOfPrism/2)*np.cos(deviationAngle/2+angleOfPrism/2)/2-np.sin(deviationAngle/2+angleOfPrism/2)*math.cos(angleOfPrism/2)/2)/(math.sin(angleOfPrism/2)))**2)

# Regressionn
DataSet = odr.RealData(meanInverseSquareWavelengths, refractiveIndex, standardDeviationInverseSquareWavelengths, uncertaintyRefractiveIndex)
ODRInstantiate = odr.ODR(DataSet, odr.unilinear)
Value = ODRInstantiate.run()

B, A = Value.beta
sdB, sdA = Value.sd_beta

#Tabulating Data
print("Wavelengths Data:")
print(tabulate(wavelengthMatrix, headers=['V1 Left','V2 Left','V1 Right','V2 Right'], showindex=colorsData))
print('Angle of prism', math.degrees(angleOfPrism), ' deg')
print("\nCalculation")
print(tabulate({
    "deviation": np.degrees(deviationAngle),
    "lamda":meanWavelengths,
    "s.d. lambda":standardDeviationWavelengths,
    "1/lamda^2": meanInverseSquareWavelengths,
    "s.d. 1/lambda^2": standardDeviationInverseSquareWavelengths,
    "mu": refractiveIndex,
    "delta mu": uncertaintyRefractiveIndex
}, showindex=colorsData, headers='keys'))

print('\nRegression Data')
print("B = ", B,"+/-", sdB)
print("A = ", A,"+/-", sdA)

#Calculating d from standard wavelengths
sineThetas = np.sin(wavelengthData)
test_d = colorsWavelength/np.mean(sineThetas, axis=1)
gratings = 0.02254/test_d

meanGratings = round(np.mean(gratings))
errorGratings = round(np.std(gratings)/math.sqrt(len(gratings)))

print("Analysis 1:")
print(tabulate({
    "Colors": colorsData,
    "True Wavelengths": colorsWavelength,
    "sin(theta)":np.mean(sineThetas, axis=1),
    "Gratings": gratings
}, headers='keys'))
print("Grating:", meanGratings,"+/-", errorGratings)


#Graphing
xFit = np.linspace(3.5e12,8e12,1000)
yFit = B*xFit + A

plt.scatter(meanInverseSquareWavelengths, refractiveIndex)
plt.errorbar(meanInverseSquareWavelengths, refractiveIndex, xerr=standardDeviationInverseSquareWavelengths, yerr=uncertaintyRefractiveIndex, linestyle='', capsize=3)
plt.plot(xFit, yFit, label='Regression curve')
plt.title(r"$\mu_{\lambda}\text{ = A + }\frac{B}{\lambda^{2}}\text{ graph}$")
plt.xlabel(r'$\frac{1}{\lambda^{2}}\text{ in }m^{-2}$')
plt.ylabel(r'$\mu_{\lambda}$')
plt.legend(loc="upper right")
plt.ylim(1.425, 1.6)
plt.show()