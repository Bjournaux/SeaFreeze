import scipy.interpolate as interp
import scipy.io as sio
import numpy as np

from mlbspline import *

# A spline giving apparent sound speed in Na2SO4 solutions
# dimensions of spline are pressure, temperature, and concentration

splineFile = 'appSS_BSpline.mat'
spd = load.loadSpline(splineFile)
velcML = sio.loadmat('velc.mat')['velc']

P = np.log10(np.logspace(-1,3.7,20))
T = np.linspace(250,800,10)
X = np.linspace(0,3,3)
x = np.array([P,T,X])
# x = np.array([-.95,285,2.1])
y = eval.evalMultivarSpline(spd,x)

maxDiff = abs((y - velcML)).max()
print('The maximum difference between the Matlab calculated spline values and those calculated by ' +
      'evalMultiVarSpline is ' + str(maxDiff))




