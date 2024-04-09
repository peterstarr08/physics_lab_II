import numpy as np
import math

#uncertainities
uMass = 1 #in gm
uTemp = 0.1 #in deg C

#Equivalent mass of flask and water
equivalentMassFlask = 23
waterMass = 183

#Aluminium
mass_Al = 44
steamChamberTemperature_Al = 94.0
waterTemp_Al = np.array([24.2, 28.5])


#Beads
mass_Beads = 58
steamChamberTemperature_Beads = 99.0
waterTemp_Beads = np.array([26.2, 29.0])

#Copper
mass_Cu = 107
steamChamberTemperature_Cu = 94.0
waterTemp_Cu = np.array([27.1, 30.3])


def analyse(massShots, steamChamberTemp, waterTemp):
    c = 1*((waterMass+equivalentMassFlask)*(waterTemp[1]-waterTemp[0]))/(massShots*(steamChamberTemp-waterTemp[1]))
    dc_dmw = (waterTemp[1]-waterTemp[0])/(steamChamberTemp-waterTemp[1])*1/massShots
    dc_dmf = dc_dmw

    dc_dtw = (waterMass+equivalentMassFlask)/(steamChamberTemp-waterTemp[1])*1/massShots
    dc_dms = -1*(waterMass+equivalentMassFlask)*(waterTemp[1]-waterTemp[0])/(steamChamberTemp-waterTemp[1])*1/(massShots**2)

    dc_dts = -1*(waterMass+equivalentMassFlask)*(waterTemp[1]-waterTemp[0])/((steamChamberTemp-waterTemp[1])**2)*1/massShots
    dc_dtm = (-1*(massShots*(steamChamberTemp-waterTemp[1])*(waterMass+equivalentMassFlask)) - (massShots*(waterMass+equivalentMassFlask)*(waterTemp[1]-waterTemp[0])))/((massShots*(steamChamberTemp-waterTemp[1]))**2)

    uncertaintyC = math.sqrt(math.pow(dc_dmw*uMass,2)+
                             math.pow(dc_dmf*uMass,2)+
                             math.pow(dc_dtw*uTemp, 2)+
                             math.pow(dc_dms*uMass, 2)+
                             math.pow(dc_dts*uTemp, 2)+
                             math.pow(dc_dtm*uMass, 2)
                             )

    print(c, '+/-', uncertaintyC, 'cal/g/degC')
print("Aluminium Shots:")
analyse(mass_Al, steamChamberTemperature_Al, waterTemp_Al)

print("\nWhite beads: ")
analyse(mass_Beads, steamChamberTemperature_Beads, waterTemp_Beads)

print("\nCopper Shots")
analyse(mass_Cu, steamChamberTemperature_Cu, waterTemp_Cu)