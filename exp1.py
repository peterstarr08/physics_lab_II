import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#Data
temperature = np.array([298, 308, 313, 318, 323, 328, 333, 338, 343, 348])
resistanceRising = np.array([10, 5.63, 4.54, 3.77, 3.184, 2.339, 2.11, 1.839, 1.571, 1.513])
resistanceFalling = np.array([10, 6.44, 5.43, 4.48, 3.627, 2.993, 2.473, 1.978, 1.667, 1.411])

#Creating metal data
sampleSize = len(temperature)

#Creating empty 1-D numpy array to store average of rising and falling resistance
resistanceAverage = np.empty(sampleSize)
for i in range(0, len(temperature)):
    resistanceAverage[i] = (resistanceFalling[i] + resistanceRising[i])/2


for i in range(0, sampleSize-1):
    R1, R2 = resistanceRising[i], resistanceRising[i+1]
    T1, T2 = temperature[i], temperature[i+1]
    
    num = math.log(R1)/math.log(R2)
    denom = 1/T1 - 1/T2

    print(num/denom)