import numpy as np
import math
import sys

rhoVec = np.arange(0.01,0.8,0.001)

oMat = []

rMin = 2**(1./6)
rc = 2.5
T = 1.

for i in range (len(rhoVec)):
	rho = rhoVec[i]
	eta = math.pi*rho

	pid = 1
	phs = (1+4*rho+4*(rho**2)-4*(rho**3)+rho**4)/((1-rho)**4)
	pdisp_long = 16*math.pi*rho*(1./T)*(-(1/(9*(rc**9)))+(1/(9*(rMin**9)))+(1/(3*(rc**3)))-(1/(3*(rMin**3))))
	pdisp_short = -(4./3)*math.pi*(rMin**3)*rho*(1./T)
	pdisp = pdisp_long+pdisp_short

	oMat.append([rho,pid+phs+pdisp])

np.savetxt('spinodal.xvg',oMat,fmt="%5f")

