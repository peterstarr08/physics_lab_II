import numpy as np
import matplotlib.pyplot as plt

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
    #R vs T scatter
    plt.scatter(temperatureData, resistanceData, s=5)
    #R vs T error bars for both resistance and temperature
    plt.errorbar(temperatureData, resistanceData, rU, tU, linestyle='')
    #formatting the plot
    plt.title(title+ r"$\text{R}(k\Omega)\text{ vs }T(K)$")
    plt.ylabel(r"$\text{Resistance}(k\Omega)$")
    plt.xlabel("Temperature(K)")
    print('\n\n', title)

    #polynomial setup with constant, degree 1 & 3 
    poly_model = odr.polynomial([1,3])
    
    #calculating error in 1/T & ln(R)
    errorTemp = []
    for val in temperatureData: 
        errorTemp.append(tU/(val*val))
    errorResistence = []
    for val in resistanceData:
        errorResistence.append(rU/(val))

    #setting up data model with uncertainty data
    dataSet = odr.RealData(np.log(resistanceData), np.reciprocal(temperatureData), errorResistence, errorTemp)
    #setting up the required data, cubic model and running it
    odrInstantiate = odr.ODR(dataSet, poly_model)
    val = odrInstantiate.run()

    #getting regression parameters with their standard deviation data
    A,B,C = val.beta
    a,b,c = val.sd_beta
    print("A=", A, '+/-', a)
    print("B=", B, '+/-', b)
    print("C=", C, '+/-', c)

    #adding regression curve
    yFit = np.linspace(0.6, 20, 1000)    
    xFit = 1/(A+B*np.log(yFit)+C*(np.log(yFit)**3))

    #more formatting and saving figure
    plt.plot(xFit, yFit, label="Regression curve")
    plt.legend(loc="upper right")
    plt.savefig(title+' R vs T.png', dpi=300)

    plt.clf()
    #Scatter for 1/T vs ln(R(T))
    plt.scatter(np.log(resistanceData),np.reciprocal(temperatureData),  s=5)
    #1/T vs ln(R(T)) error bars
    plt.errorbar(np.log(resistanceData),np.reciprocal(temperatureData), errorTemp, errorResistence, linestyle='')
    #Plot formatting
    plt.title(title+" 1/T vs ln(R)")
    plt.ylabel("1/T")
    plt.xlabel("ln(R)")

    #adding regression curve
    xFit = np.linspace(-0.2, 2.5, 1000)
    yFit = A+B*(xFit)+C*(xFit**3)
    plt.plot(xFit, yFit, label="Regression curve")

    #formatting and savinf figure
    plt.legend(loc="upper right")
    plt.savefig(title+' 1_T vs ln(R).png', dpi=300)


analysis(celciusData, 1, risingTemperatureResistance, 0.001, "Rising Temperature Data- ")
analysis(celciusData, 1, loweringTemperatureResistance, 0.001, "Cooling Temperature Data- ")
