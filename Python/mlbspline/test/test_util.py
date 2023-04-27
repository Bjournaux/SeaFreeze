from numpy import finfo, where

eps = finfo(float).eps

def compareToMatlabOutput(mlbsplineOutput, matlabOutput, measureName):
    absDiff = abs(mlbsplineOutput - matlabOutput)
    relDiff = absDiff / where(matlabOutput != 0, matlabOutput, eps)
    print("\n".join(['Max difference between Matlab and mlbspline out for ' + measureName + ':',
          '\t\tabsolute:\t\t' + str(absDiff.max()), '\t\trelative:\t\t'+ str(relDiff.max())]))
