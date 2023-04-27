from scipy.io import loadmat
import numpy as np

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# from pylab import get_current_fig_manager

from mlbspline import *

# A spline giving Gibbs energy in water
# dimensions of spline are P,T

splineFile = 'sp_IAPWS.mat'
spd = load.loadSpline(splineFile,'sp_IAPWS')
gibbsML = loadmat('IAPWS_Gibbs.mat')['G']

P = np.arange(.1,1000,1) #MPa
T = np.arange(240,501) #K
x = np.array([P,T])
y = eval.evalMultivarSpline(spd,x)

maxDiff = abs((y - gibbsML)).max()
print('The maximum difference between the Matlab calculated spline values and those calculated by ' +
      'eval is ' + str(maxDiff))


# fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.set_xlabel('Pressure (MPa)',labelpad=20)
# ax.set_ylabel('Temperature (K)',labelpad=20)
# ax.set_zlabel('Gibbs Energy (J/Mg)',labelpad=20)
# P, T = np.meshgrid(P,T)
# ax.plot_surface(P,T,y.T)
# plt.title('IAPWS Gibbs Energy for Pure Water')
# plt.show()

# check derivatives
gibbsP1 = loadmat('IAPWS_Gibbs.mat')['Gd1P']
yp1 = eval.evalMultivarSpline(spd,x,[1, 0])
maxDiff = abs((yp1 - gibbsP1)).max()
print('For the first deriv wrt P, maximum difference between the Matlab calculated spline values and those calculated by ' +
      'eval is ' + str(maxDiff))


gibbsT1 = loadmat('IAPWS_Gibbs.mat')['Gd1T']
yt1 = eval.evalMultivarSpline(spd,x,[0, 1])
maxDiff = abs((yt1 - gibbsT1)).max()
print('For the first deriv wrt T, maximum difference between the Matlab calculated spline values and those calculated by ' +
      'eval is ' + str(maxDiff))

gibbsP2 = loadmat('IAPWS_Gibbs.mat')['Gd2P']
yp2 = eval.evalMultivarSpline(spd,x,[2, 0])
maxDiff = abs((yp2 - gibbsP2)).max()
print('For the second deriv wrt P, maximum difference between the Matlab calculated spline values and those calculated by ' +
      'eval is ' + str(maxDiff))

gibbsT2 = loadmat('IAPWS_Gibbs.mat')['Gd2T']
yt2 = eval.evalMultivarSpline(spd,x,[0, 2])
maxDiff = abs((yt2 - gibbsT2)).max()
print('For the second deriv wrt T, maximum difference between the Matlab calculated spline values and those calculated by ' +
      'eval is ' + str(maxDiff))

gibbsPT = loadmat('IAPWS_Gibbs.mat')['GdPT']
ypt = eval.evalMultivarSpline(spd,x,[1, 1])
maxDiff = abs((ypt - gibbsPT)).max()
print('For the first derivative wrt both P and T, maximum difference between the Matlab calculated spline values and those calculated by ' +
      'eval is ' + str(maxDiff))





