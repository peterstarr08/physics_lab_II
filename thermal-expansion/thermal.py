import numpy as np
import math
from scipy import odr
import matplotlib.pyplot as plt

DIAL_LC = 0.01
TEMP_LC = 0.1
LENGTH_LC = 1

#61-63

#Brass
Brass_Length = 70.1*10
Brass_Readings = np.array([
    [42.9, 20],
    [42.3, 19],
    [41.7, 18],
    [40.8, 17],
    [40.3, 16],
    [39.7, 15],
    [38.9, 14],
    [38.0, 13],
    [37.0, 11],
    [36.1, 10],
    [35.3, 9],
    [34.6, 8],
    [33.6, 7],
    [32.9, 6],
    [32.1, 5],
    [31.2, 4],
    [30.3, 3],
    [29.6, 2],
    [28.8, 1],
    [27.7, 0]
])
Brass_Readings[:,0] = Brass_Readings[:,0]*1.2
#Copper
Copper_Length = 70.3*10
Copper_Readings = np.array( [
    [56.7, 40],
    [56.0, 39],
    [55.3, 38],
    [54.3, 37],
    [53.6, 36],
    [52.6, 35],
    [52.3, 34],
    [51.7, 33],
    [50.5, 32],
    [49.3, 29],
    [48.5, 28],
    [48.2, 27],
    [47.9, 26],
    [47.5, 25],
    [46.2, 24],
    [45.8, 23],
    [44.9, 22],
    [44.2, 21],
    [43.6, 20],
    [43.1, 19],
    [42.6, 18],
    [41.8, 17],
    [41.1, 16],
    [40.4, 15],
    [39.7, 14],
    [38.9, 13],
    [38.1, 12],
    [37.6, 11],
    [36.8, 10],
    [36.2, 9],
    [35.5, 8],
    [34.5, 7],
    [33.9, 6],
    [32.9, 5],
    [32.3, 4],
    [31.5, 3],
    [30.8, 2],
    [29.9, 1],
    [29.1, 0]
])
Copper_Readings[:,0]=Copper_Readings[:,0]*1.2
#Aluminium
Aluminium_Length = 70.4*10
Aluminium_Readings = np.array([
    [60.9, 45],
    [60.4, 44],
    [59.1, 43],
    [58.8, 42],
    [58.3, 41],
    [57.3, 40],
    [56.8, 39],
    [55.8, 38],
    [55.3, 37],
    [54.4, 36],
    [54.0, 35],
    [53.1, 34],
    [52.4, 33],
    [51.7, 32],
    [50.6, 31],
    [50.3, 30],
    [49.5, 29],
    [48.9, 28],
    [47.8, 27],
    [47.2, 26],
    [46.5, 25],
    [45.9, 24],
    [44.9, 23],
    [44.1, 22],
    [43.5, 21],
    [42.9, 20],
    [41.9, 19],
    [41.2, 18],
    [40.8, 17],
    [39.9, 16],
    [39.2, 15],
    [38.6, 14],
    [37.9, 13],
    [37.1, 12],
    [36.4, 11],
    [35.6, 10],
    [34.8, 9],
    [34.2, 8],
    [33.7, 7],
    [32.9, 6],
    [32.2, 5],
    [31.5, 4],
    [30.8, 3],
    [30.3, 2],
    [29.4, 1],
    [28.9, 0]
])
Aluminium_Readings[:,0] = Aluminium_Readings[:,0]*1.2

#Analysis function
def analysis(initialLength ,dialData, title):
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
    plt.title(title+' - '+r'$\Delta\text{L = }\alpha\text{L}_{0}T\text{ - }\alpha\text{L}_{0}T_{0} $')
    plt.ylabel(r'$\Delta L\text{( in mm)}$')
    plt.xlabel(r'$\text{Temperature T}(\degree C)$')
    plt.plot(xFit, yFit, label='Regression curve')
    plt.show()

analysis(Brass_Length, Brass_Readings,'Brass')
analysis(Aluminium_Length, Aluminium_Readings, 'Steel')
analysis(Copper_Length, Copper_Readings, 'Copper')