import numpy as np
import matplotlib.pyplot as plt
import math

from scipy import odr

#Celcius data starting from 23, 24, 26, 28,....88
celciusData = np.insert(np.array(range(297, 363, 2), dtype=float), 0, 296)

#Resistance in kilo ohms whe temperature rised corresponding the celcius data
risingTemperatureResistance = np.array([11.55, 10.60, 9.63, 8.73, 8.05, 7.28, 6.64, 6.09, 5.59,
                                        5.13, 4.72, 4.32, 3.95, 3.64, 3.350, 3.064, 2.799, 2.579,
                                        2.379, 2.214, 2.044, 1.874, 1.745, 1.619, 1.515, 1.408,
                                        1.329, 1.253, 1.180, 1.110, 1.048, 0.999, 0.949, 0.897])
loweringTemperatureResistance = np.array([10.72, 10.45, 9.74, 9.11, 8.60, 8.16, 7.77, 7.25, 6.64,
                                          6.07, 5.61, 5.07, 4.45, 3.968, 3.498, 3.217, 2.942,
                                           2.734, 2.508, 2.332, 2.168, 2.008, 1.876, 1.720,
                                           1.646, 1.553, 1.458, 1.356, 1.253, 1.159, 1.153, 1.009,
                                           0.973, 0.917])

def analysis(temperatureData, tU, resistanceData, rU, title):
    plt.clf()
    plt.scatter(temperatureData, resistanceData)
    plt.errorbar(temperatureData, resistanceData, tU, rU, linestyle='')
    plt.title(title+ r"$\text{R}(k\Omega)\text{ vs }T(K)$")
    plt.ylabel(r"$\text{Resistance}(k\Omega)$")
    plt.xlabel("Temperature(K)")
    print('\n\n', title)
    poly_model = odr.polynomial([1,3])
    dataSet = odr.RealData(np.log(resistanceData), np.reciprocal(temperatureData), rU, tU)
    
    odrInstantiate = odr.ODR(dataSet, poly_model)
    val = odrInstantiate.run()
    A,B,C = val.beta
    a,b,c = val.sd_beta
    print("A=", A, '+/-', a)
    print("B=", B, '+/-', b)
    print("C=", C, '+/-', c)

    yFit = np.linspace(0.5, 20, 1000)
    
    xFit = 1/(A+B*np.log(yFit)+C*(np.log(yFit)**3))

    plt.plot(xFit, yFit, label="Regression curve")
    plt.legend(loc="upper right")

    plt.savefig(title+'.png', dpi=300)
    # plt.show()


analysis(celciusData, 1, risingTemperatureResistance, 0.001, "Rising Temperature Data- ")
analysis(celciusData, 1, loweringTemperatureResistance, 0.001, "Cooling Temperature Data- ")
