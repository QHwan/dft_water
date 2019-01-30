import numpy as np
import sys
import math

rhoVec = np.arange(0,1.001,0.001)

g1Vec = []
g2Vec = []

for i in range (len(rhoVec)):
	rho = rhoVec[i]

	g = 1/(1-math.pi*rho/6.) + (1/4.)*((math.pi*rho)/(pow(1-math.pi*rho/6,2.))) + (1/72.)*(pow(math.pi*rho,2.))/(pow(1-math.pi*rho/6,3.))
	gp = (5*math.pi/12.)/(pow(1-math.pi*rho/6.,2.)) + (math.pi*math.pi/9.)*(rho)/(pow(1-math.pi*rho/6,3.)) + (math.pi*math.pi*math.pi/144.)*(pow(rho,2.))/(pow(1-math.pi*rho/6,4.))

	

	g1Vec.append(g)
	g2Vec.append(gp)

g1Vec = np.gradient(g1Vec,0.001)


o1Vec = []
o1Vec.append(rhoVec); o1Vec.append(g1Vec); o1Vec.append(g2Vec)
o1Vec = np.transpose(o1Vec)




np.savetxt('1.xvg',o1Vec,fmt='%5f')
