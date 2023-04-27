from scipy.io import loadmat
import numpy as np

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import get_current_fig_manager

from mlbspline import *

# A spline giving volume of pure water
# dimensions of spline are pressure (GPa) and temperature (K)
# ice VII Mie-Grüneisen equation of state based on X-Ray diffraction data from Bezacier et al. (2014)

splineFile = 'iceVII_EOS.mat'
spd = load.loadSpline(splineFile)
volML = loadmat('iceVIIvol.mat')['vol']

P = np.logspace(0,np.log10(20),50)
T = np.linspace(0,550,100)
x = np.array([P,T])
y = eval.evalMultivarSpline(spd,x)

maxDiff = abs((y - volML)).max()
print('The maximum difference between the Matlab calculated spline values and those calculated by ' +
      'evalMultiVarSpline is ' + str(maxDiff))

# load data
dat = loadmat('data_iceVII.mat')['data']
# data provided is T, P, V - rearrange to match spline
dat = np.array([[d[1],d[0],d[2]] for d in dat if 0 < d[1] < 20 and 0 < d[0] < 550])



fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlabel('Pressure (GPa)',labelpad=20)
ax.set_ylabel('Temperature (K)',labelpad=20)
ax.set_zlabel('Volume ($m^3$/kg)',labelpad=20)
P, T = np.meshgrid(P,T)
ax.plot_surface(P,T,y.T)
ax.scatter(dat[:,0],dat[:,1],dat[:,2],c='k')
plt.title('Ice VII Volume by Pressure and Temperature')
plt.show()

