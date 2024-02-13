import numpy as np
import math
from scipy import odr
import matplotlib.pyplot as plt

DIAL_LC = 0.01
TEMP_LC = 0.1
LENGTH_LC = 1

#Brass
Brass_Length = 57.8*10
Brass_Readings = np.array([
    [30.1, 7],
    [28.8, 6],
    [27.7, 5],
    [26.7, 4],
    [25.8, 3],
    [25, 2],
    [23.8, 1],
    [22.8, 0]
])

#Copper
Copper_Length = 57.9*10
Copper_Readings = np.array([
    [78.8, 72],
    []
])

#Aluminium
Aluminium_Length = 57.7*10
Aluminium_Readings = np.array([
    [77.9, 73],
    []
])

#Analysis function
def analysis(initialLength ,dialData):
    temperatureData = dialData[:,0]
    lengthData = dialData[:,1]*DIAL_LC

    #Calculation
    DataSet = odr.RealData(temperatureData, lengthData, TEMP_LC, DIAL_LC)
    ODRInstantiate = odr.ODR(DataSet, odr.unilinear)
    Value = ODRInstantiate.run()

    m,c  = Value.beta
    _m, _c = Value.sd_beta
    
    alpha = m/initialLength
    _alpha = math.sqrt((_m/initialLength)**2+(m*LENGTH_LC/(initialLength**2))**2)

    print('Coefficient of expansion: ', alpha,'+/-',_alpha, '/degC')

    yFit = np.linspace((c-1)*DIAL_LC, (dialData[0,1]+1)*DIAL_LC, 1000)
    xFit = (yFit-c)/m

    plt.scatter(temperatureData, lengthData, s=5)
    plt.errorbar(temperatureData, lengthData, xerr=TEMP_LC, yerr=DIAL_LC, linestyle='')
    plt.plot(xFit, yFit)
    plt.show()

analysis(Brass_Length, Brass_Readings)