import math
import numpy as np
from scipy.stats import tstd


#gravity
g = 9.8

#distance between parallel plates
d = 5e-3

#distance between the chosen top and bottom lines on the monitor
L = 1e-3

#desnity of the oil
rho_oil = 928

#density of the air
rho_air  = 1

#room temperature
roomTemp = 22 + 273

#atmospheric pressure(in m of Hg)
P = 0.76

#coefficient of viscosity of air
eta = 1.8333e-5

#correct factor c
c = 6.17e-8

#C calculation
C = 4*math.pi*d*g*(rho_oil - rho_air)/3
print('C', C)

#D calculation
D = 9*eta/(2*g*(rho_oil-rho_air))
print("D", D)

#zeta
zeta = c/(2*P)
print("zeta", zeta)

#Dynamic data
dynamicDrop = np.array([
    [np.mean([7.98, 8.04, 7.88, 7.75]), np.mean([6.67, 7.21, 7.88, 7.1, 7.34]), 120],
    [np.mean([13.95, 13.39, 13.62]), np.mean([10.68, 11.23, 11.14]), 120.0],
    [np.mean([8.59, 8.86]), np.mean([13.3, 13.37]), 240.0],
    [np.mean([12.64, 12.07, 12.19]), np.mean([3.66, 3.73, 3.77]), 240.0],
    [np.mean([4.64, 4.59, 4.51, 4.38]), np.mean([6.02, 6.03, 6.09, 6.05]), 170.0],
    [np.mean([6.27, 6.36, 6.33]), np.mean([5.3, 5.65, 5.63]), 170.0],
    [np.mean([8.94, 8.25]), np.mean([16.84, 18.11]), 120.0]
])

v_f_d = L/dynamicDrop[:,0]

xi_d = D*v_f_d
r_d = -zeta + np.sqrt(zeta*zeta + xi_d)
r3_d = np.power(r_d, 3)
T_d = 1 + dynamicDrop[:, 0]/ dynamicDrop[:, 1]
ne = C*np.divide(np.multiply(T_d,r3_d), dynamicDrop[:,2])

min_ne_d = np.min(ne)
div_ne_d = ne/min_ne_d
neff_d = np.round(div_ne_d)
ne_neff_d = ne/neff_d 
print('\n\nDynamic Data',dynamicDrop)
print('Mean velcoity', v_f_d)
print('xi',xi_d)
print('r', r_d)
print('r3', r3_d)
print('T', T_d)
print('ne', ne)

print('ne divided by the lowest', div_ne_d)
print('n_eff', neff_d)
print('ne/n_eff', ne_neff_d)
print("e=", np.mean(ne_neff_d), 's. error=',tstd(ne_neff_d)/(len(ne_neff_d)**0.5))

#Balancing method
balanceDrop = np.array([
    [np.mean([9.64, 9.34]), np.mean([165, 166])],
    [np.mean([21.75, 21.4]), np.mean([120, 120])],
    [np.mean([18.19, 18.19]), np.mean([124, 124])],
    [np.mean([17.31, 17.61]), np.mean([151, 149])],
    [np.mean([22.03, 21.43]), np.mean([201, 205])],
    [np.mean([21.71, 19.87]), np.mean([196, 200])],
    [np.mean([16.1, 15.81]), np.mean([165, 160])],
])

v_f_b = L/balanceDrop[:, 0]
xi_b = D*v_f_b
r_b = -zeta + np.sqrt(zeta*zeta + xi_b)
r3_b = np.power(r_b, 3)
ne_b = C*1*r3_b/balanceDrop[:, 1]
min_ne_b = np.min(ne_b)
div_ne_b = ne_b/min_ne_b
neff_b = np.round(div_ne_b)
ne_neff_b = ne_b/neff_b

print('\n\nBalance Drop', balanceDrop)
print('Mean Velocity', v_f_b)
print('xi',xi_b)
print('r', r_b)
print('r3', r3_b)
print('ne', ne_b)

print('ne divided by the lowest', div_ne_b)
print('n_eff', neff_b)
print('ne/n_eff', ne_neff_b)
print('e=', np.mean(ne_neff_b), 's. error=', tstd(ne_neff_b)/(len(ne_neff_b)*0.5))
