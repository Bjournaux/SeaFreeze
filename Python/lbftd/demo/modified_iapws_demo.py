import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D     # even though this is not used directly, it is required to do 3d plots
from lbftd import statevars, loadGibbs as lg, evalGibbs as eg

# load and evaluate the spline
# this is a high-pressure spline that starts with IAPWS 95 below 1 GPa
# but modifies it above 1 GPa to better match DAC and shockwave data
water_spline = lg.loadGibbsSpline('water_demo_iapws_modified.mat')

P = np.linspace(0, 2000, num=200)       # pressure, in MPa
T = np.linspace(200, 1000, num=200)        # temperature, in K
# requested thermodynamic state variables
# (see README for full list of available state vars or don't list any to get the full set):
# - rho: density in kg m^-3
# - Cp: isobaric specific heat in J kg^-1 K^-1
# - Kt: isothermal bulk modulus in MPa
# - alpha:  thermal expansivity in K-1
tdstate = eg.evalSolutionGibbsGrid(water_spline['sp'], np.array([P, T]), 'rho', 'Cp', 'Kt', 'alpha', failOnExtrapolate=False)
# for full set of implemented statevars: eg.evalSolutionGibbsGrid(water_spline['sp'], np.array([P, T]))

# graph the output
rcParams['axes.labelpad'] = 10 # add some padding to prevents axis labels from covering ticks
fig = plt.figure()
pP, pT = np.meshgrid(P, T)

rho_ax = fig.add_subplot(221, projection='3d')
rho_ax.set_xlabel('Pressure ($MPa$)')
rho_ax.set_ylabel('Temperature ($K$)')
rho_ax.set_zlabel('Density ($kg\: m^{-3}$)')  # use LaTex coding in labels
rho_surf = rho_ax.plot_surface(pP, pT, tdstate.rho)
rho_ax.invert_yaxis()
rho_ax.set_zlim3d(-5000, 10000)

cp_ax = fig.add_subplot(222, projection='3d')
cp_ax.set_xlabel('Pressure ($MPa$)')
cp_ax.set_ylabel('Temperature ($K$)')
cp_ax.set_zlabel('Specific Heat ($J\: kg^{-1}\: K^{-1}$)')
cp_surf = cp_ax.plot_surface(pP, pT, tdstate.Cp)
cp_ax.invert_yaxis()

kt_ax = fig.add_subplot(223, projection='3d')
kt_ax.set_xlabel('Pressure ($MPa$)')
kt_ax.set_ylabel('Temperature ($K$)')
kt_ax.set_zlabel('Isothermal Bulk Modulus ($MPa$)')
kt_surf = kt_ax.plot_surface(pP, pT, tdstate.Kt)
kt_ax.invert_yaxis()
# IAPWS has spike for this state var.  Limit the range of z-values shown
kt_ax.set_zlim3d(-3000, 10000)

alpha_ax = fig.add_subplot(224, projection='3d')
alpha_ax.set_xlabel('Pressure ($MPa$)')
alpha_ax.set_ylabel('Temperature ($K$)')
alpha_ax.set_zlabel('Thermal Expansivity ($K^{-1}$)')
alpha_surf = alpha_ax.plot_surface(pP, pT, tdstate.alpha)
alpha_ax.invert_yaxis()

fig.suptitle('IAPWS 95 below 1 GPa \n modified above 1 GPa to match DAC and shockwave data')

plt.show()









