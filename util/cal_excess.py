import numpy as np
import sys
import math

iFileName = sys.argv[1]

iMat = np.loadtxt(iFileName)
iMat = np.transpose(iMat)

rVec = iMat[0]
dVec = iMat[1]

dv = dVec[len(dVec)-1]

for i in range (len(dVec)):
	dVec[i] = dVec[i] - dv
	dVec[i] = dVec[i] * 4*math.pi*rVec[i]*rVec[i]

nexc = np.trapz(dVec,rVec)

print nexc

