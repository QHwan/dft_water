import numpy as np

zVec = np.arange(0.1,0.5,0.001)
fVec = np.zeros(len(zVec))
gVec = np.zeros(len(zVec))

for i in range (len(zVec)):
	z = zVec[i]
	fVec[i] += z*(8-9*z+3*z*z)/((1-z)**3) + 1/z
	gVec[i] += (4-z)/((1-z)**4)
	
f2Vec = np.gradient(fVec,0.001)

oMat = []
oMat.append(zVec)
oMat.append(f2Vec)
oMat.append(gVec)

np.savetxt('text.xvg',np.transpose(oMat),fmt='%5f')



