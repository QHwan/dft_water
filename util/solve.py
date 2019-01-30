import numpy as np
import sys
import math
def fp(eta,T,rMin,rc):
	rho = 6*eta/PI
	p = (rho)*((1+eta+eta*eta-eta*eta*eta)/(math.pow((1-eta),3))) \
	    + (48*eta*rho)/(T)*(-1/(9*math.pow(rc,9))+1/(9*math.pow(rMin,9))+1/(3*math.pow(rc,3))-1/(3*math.pow(rMin,3))) \
	    - (4*rho*eta/T)*(math.pow(rMin,3))
	return p

def fmu(eta,T,rMin,rc):
	rho = 6*eta/PI
	#mu = math.log(eta*math.pow(1/T,1.5)) + (eta*(8-9*eta+3*eta*eta))/(math.pow(1-eta,3)) + (96*eta)/(T)*(-1/(9*math.pow(rc,9))+1/(9*math.pow(rMin,9))+1/(3*math.pow(rc,3))-1/(3*math.pow(rMin,3))) - 8/T*eta*(math.pow(rMin,3))
	mu = math.log(eta) \
	     + (eta*(8-9*eta+3*eta*eta))/(math.pow(1-eta,3)) \
	     + (96*eta)/(T)*(-1/(9*math.pow(rc,9))+1/(9*math.pow(rMin,9))+1/(3*math.pow(rc,3))-1/(3*math.pow(rMin,3))) \
	     - (8*eta/T)*(math.pow(rMin,3))

	return mu


oFileName = 'phase_wca.xvg'
PI = 3.1415
TVec = np.arange(1.0,1.0+0.01,0.2)
denVec = np.arange(0.001,1.0,0.0005)
rc = 2.5
rMin = 2.**(1/6.)
T = 1
tol = 1E-3

pMat = []
oMat = []

			


#for den in denVec:
#	pMat.append([den,fp(den*PI/6,T,rMin,rc)])
#	oMat.append([den,fmu(den*PI/6,T,rMin,rc)])
#np.savetxt("p.xvg",pMat,fmt='%5f')
#np.savetxt("mu.xvg",oMat,fmt='%5f')
#exit(1)




for T in TVec:
	pMat = []
	for i in range (len(denVec)-1):
		d1 = denVec[i]
		p1 = fp(PI/6*d1,T,rMin,rc)

		d2 = denVec[i+1]
		p2 = fp(PI/6*d2,T,rMin,rc)

		if p2-p1 > 0: # stable condition
			for j in range (len(denVec)-1, i+100, -1):
				d3 = denVec[j]
				p3 = fp(PI/6*d3,T,rMin,rc)

				d4 = denVec[j-1]
				p4 = fp(PI/6*d4,T,rMin,rc)

				if p3 - p4 > 0:
					if abs(p3 - p1) < tol:
						print p3-p1
						pMat.append([d1,d3])
				else:
					break
		else:
			pass

	diffVec = []
	for i in range (len(pMat)-1):
		d1 = pMat[i][0]
		d2 = pMat[i][1]
		mu1 = fmu(PI/6*d1,T,rMin,rc)
		mu2 = fmu(PI/6*d2,T,rMin,rc)

		diffVec.append(abs(mu2-mu1))

	index = np.argmin(diffVec)
	print diffVec[index]

	mu1 = fmu(pMat[index][0]*PI/6,T,rMin,rc)
	mu2 = fmu(pMat[index][1]*PI/6,T,rMin,rc)
	mu3 = (mu1+mu2)/2

#print pMat[index][0], pMat[index][1], diffVec[index]
	oMat.append([pMat[index][0],pMat[index][1],mu1,mu2,(mu1+mu2)/2,T])
			

np.savetxt(oFileName,oMat,fmt='%5f')


