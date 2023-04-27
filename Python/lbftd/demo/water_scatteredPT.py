from random import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D
from lbftd import statevars, loadGibbs as lg, evalGibbs as eg


# generate n scattered data points
n = 500
minT = 250;     maxT = 350
minP = 0;       maxP = 1000
toP = lambda r: (maxP - minP)*r + minP
toT = lambda r: (maxT - minT)*r + minT
PT = np.empty(n, np.object)
for i in np.arange(0, n):
    PT[i] = (toP(random()), toT(random()))

water_spline = lg.loadGibbsSpline('water_demo_spline.mat')
tdstate = eg.evalSolutionGibbsScatter(water_spline['sp'], PT, 'rho', 'Cp', 'Kt', 'alpha')

# graph the output
rcParams['axes.labelpad'] = 10 # add some padding to prevents axis labels from covering ticks
fig = plt.figure()
# sorted sets of P and T for creating the meshgrid for plotting
P = [pt[0] for pt in PT]
T = [pt[1] for pt in PT]


rho_ax = fig.add_subplot(221, projection='3d')
rho_ax.set_xlabel('Pressure ($MPa$)')
rho_ax.set_ylabel('Temperature ($K$)')
rho_ax.set_zlabel('Density ($kg\: m^{-3}$)')  # use LaTex coding in labels
rho_surf = rho_ax.scatter(P, T, tdstate.rho)
rho_ax.invert_yaxis()

cp_ax = fig.add_subplot(222, projection='3d')
cp_ax.set_xlabel('Pressure ($MPa$)')
cp_ax.set_ylabel('Temperature ($K$)')
cp_ax.set_zlabel('Specific Heat ($J\: kg^{-1}\: K^{-1}$)')
cp_surf = cp_ax.scatter(P, T, tdstate.Cp)
cp_ax.invert_yaxis()

kt_ax = fig.add_subplot(223, projection='3d')
kt_ax.set_xlabel('Pressure ($MPa$)')
kt_ax.set_ylabel('Temperature ($K$)')
kt_ax.set_zlabel('Isothermal Bulk Modulus ($MPa$)')
kt_surf = kt_ax.scatter(P, T, tdstate.Kt)
kt_ax.invert_yaxis()

alpha_ax = fig.add_subplot(224, projection='3d')
alpha_ax.set_xlabel('Pressure ($MPa$)')
alpha_ax.set_ylabel('Temperature ($K$)')
alpha_ax.set_zlabel('Thermal Expansivity ($K^{-1}$)')
alpha_surf = alpha_ax.scatter(P, T, tdstate.alpha)
alpha_ax.invert_yaxis()

plt.show()









