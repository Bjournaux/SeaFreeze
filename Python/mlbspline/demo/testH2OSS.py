import scipy.io as sio
import scipy.interpolate as interp
import numpy as np

from mlbspline import *

# A spline calculating sound speed in pure water
# dimensions of the spline are pressure and temperature

splineFile = 'sp_H2O_SS.mat'
wspd = load.loadSpline(splineFile)
waterML = sio.loadmat('water.mat')['water']

P = np.logspace(-1,3.7,50)
T = np.linspace(250,800,20)
x = np.array([P,T])
#waterSP = interp.bisplev(P,T,wspd)
w = eval.evalMultivarSpline(wspd,x)


maxDiff = abs((w - waterML)).max()
print('The maximum difference between the Matlab calculated spline values and those calculated by ' +
      'evalMultiVarSpline is ' + str(maxDiff))

