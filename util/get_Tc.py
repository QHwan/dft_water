import numpy as np
import scipy.optimize
import sys
import math
import random

def get_a(T,rc):
	rMin = 2**(1./6.)
	a = 2*(math.pi*pow(rMin,3))/(3*T)
	a += (8.)*math.pi*(1/T)*(1./(9*pow(rc,9.))-1./(3*pow(rc,3)) - 1/(9*pow(rMin,9.)) +1/(3*pow(rMin,3)))
	return a


def fmu_id(rho,M):
	rhoM = rho/M
	return math.log(rhoM)

def fp_id(rho,M):
	rhoM = rho/M
	return rhoM

def fmu_hs(rho):
	eta = (rho)*math.pi/6
	return eta*(8-9*eta+3*eta*eta)/(pow(1-eta,3.))

def fp_hs(rho):
	eta = rho*math.pi/6
	return rho*(4*eta-2*eta*eta)/(pow(1-eta,3.))

def fmu_disp(rho,a):
	return -2*a*(rho)

def fp_disp(rho,a):
	return -a*rho*rho

################################### ass ######################################
def get_ghs(rho):
	return 1/(1-math.pi*rho/6.) + (1/4.)*((math.pi*rho)/(pow(1-math.pi*rho/6,2.))) + (1/72.)*(pow(math.pi*rho,2.))/(pow(1-math.pi*rho/6,3.))

def get_gphs(rho):
	return (5*math.pi/12.)/(pow(1-math.pi*rho/6.,2.)) + (math.pi*math.pi/9.)*(rho)/(pow(1-math.pi*rho/6,3.)) + (math.pi*math.pi*math.pi/144.)*(pow(rho,2.))/(pow(1-math.pi*rho/6,4.))

def get_D(rho,K,T_assoc):
	ghs = get_ghs(rho)
	return 4*math.pi*K*ghs*(math.exp(1/T_assoc)-1)

def get_Dp(rho,K,T_assoc):
	gphs = get_gphs(rho)
	return 4*math.pi*K*gphs*(math.exp(1/T_assoc)-1)

def get_X(rho,Ma,K,T_assoc):
	D = get_D(rho,K,T_assoc)
#	return 1/(1+math.sqrt(1+4*Ma*rho*D))
	return (-1+(1+4*Ma*rho*D)**0.5)/(2*Ma*rho*D)

def get_Xp(rho,Ma,K,T_assoc):
	D = get_D(rho,K,T_assoc)
	Dp = get_Dp(rho,K,T_assoc)
#	return -(2*Ma)/(math.sqrt(1+4*Ma*rho*D)*pow(1+math.sqrt(1+4*Ma*rho*D),2.)) * (D + rho*Dp) 
	#return 0.5*((1+4*Ma*rho*D)**-0.5)*4*Ma*Dp*(1/(2*Ma*rho*D)) + (-1+(1+4*Ma*rho*D)**0.5)/(2*Ma)*(-Dp)/(rho*rho*D*D)
	return (1/(rho*D*((1+4*Ma*rho*D)**0.5)))*(D+rho*Dp) + ((1-(1+4*Ma*rho*D)**0.5)/(2*Ma*rho*rho*D*D))*(D+rho*Dp)


def fmu_assoc(rho,Ma,K,T_assoc):
	X = get_X(rho,Ma,K,T_assoc)
	Xp = get_Xp(rho,Ma,K,T_assoc)
	if Ma == 0:
		return 0
	else:
		return Ma*(math.log(X)-0.5*X+0.5) + Ma*rho*Xp*((1/X) - 0.5)

def fp_assoc(rho,Ma,K,T_assoc):
	X = get_X(rho,Ma,K,T_assoc)
	Xp = get_Xp(rho,Ma,K,T_assoc)
	if Ma == 0:
		return 0
	else:
		return  -Ma*rho*(math.log(X) - 0.5*X + 0.5) + rho*(Ma*(math.log(X)-0.5*X+0.5) + Ma*rho*Xp*((1/X) - 0.5))
##############################################################################

def fmu_chain(rho,M):
	ghs = get_ghs(rho)
	gphs = get_gphs(rho)

	if M == 1:
		return 0
	else:
		return ((1-M)/M)*(math.log(ghs)+(rho/ghs)*(gphs))

def fp_chain(rho,M):
	ghs = get_ghs(rho)
	gphs = get_gphs(rho)

	if M == 1:
		return 0
	else:
		return  -((1-M)/M)*rho*math.log(ghs) + rho*(((1-M)/M)*(math.log(ghs)+(rho/ghs)*(gphs))
)
##############################################################################

def get_mu(rho,a,Ma,K,T_assoc,M):
	mu = fmu_id(rho,M) + M*(fmu_hs(rho) + fmu_disp(rho,a) + fmu_assoc(rho,Ma,K,T_assoc) + fmu_chain(rho,M))
	return mu


def find_peaks_num_mu(T,a,Ma,K,T_assoc,M):
	denVec = np.arange(0.001,1.5,0.001)

	test_solVec = []
	for i in range (1,len(denVec)-1):
		f = get_mu(denVec[i],a,Ma,K,T_assoc,M)
		f_before = get_mu(denVec[i-1],a,Ma,K,T_assoc,M)
		f_after = get_mu(denVec[i+1],a,Ma,K,T_assoc,M)
		if (f-f_before)*(f_after-f)<0:
			test_solVec.append(denVec[i])
	
	return len(test_solVec), test_solVec


# find mu_coex and p_coex
ofilename = sys.argv[1]

rc = 2.5

M = 1
M_vec = [1.]
Ma_vec = [4.]
epsVec = [1.,2.,3.,4.,5.,6.]
K = 1.4849*1e-4


oMat = []
for i in range (len(epsVec)):

	eps = epsVec[i]
	T_assoc = 1/eps
	Ma = 4.

	tol_ref = 1e-5
	tol = 1.
	T1 = 0.5
	T2 = 100.


	T_assoc_1 = T_assoc
	T_assoc_2 = T_assoc


	a1 = get_a(T1,rc)
	a2 = get_a(T2,rc)

	num_sol1, testVec1 = find_peaks_num_mu(T1,a1,Ma,K,T_assoc_1,M)
	num_sol2, testVec2 = find_peaks_num_mu(T2,a2,Ma,K,T_assoc_2,M)


	while (tol>tol_ref):
		T3 = (T1+T2)*0.5
		T_assoc_3 = T_assoc*T3

		a3 = get_a(T3,rc)

		num_sol3, testVec3 = find_peaks_num_mu(T3,a3,Ma,K,T_assoc_3,M)

		if num_sol1 == 2 and num_sol3 == 0:
			T2 = T3
			num_sol2 = num_sol3
			testVec2 = testVec3
		elif num_sol3 == 2 and num_sol2 == 0:
			T1 = T3
			num_sol1 = num_sol3
			testVec1 = testVec3
		else:
			print ("something is wrong")

		tol = abs(T2-T1)
		print T1, T3, T2, num_sol1, num_sol3, num_sol2

	print np.average(testVec1)

	oMat.append([1/T_assoc,(T1+T2)*0.5,np.average(testVec1)])

np.savetxt(ofilename,oMat,fmt='%5f')

