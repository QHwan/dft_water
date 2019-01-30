import numpy as np
import scipy.optimize
import sys
import math
import random

def get_a(T,rc):
	rMin = 2**(1./6.)
	a = 2*(math.pi*pow(rMin,3))/(3*T)
	a += (8.)*math.pi*(1/T)*((1/(9*pow(rc,9.))-1/(3*pow(rc,3)))-(1/(9*pow(rMin,9.))-1/(3*pow(rMin,3))))
	return a

def get_m(T,rc):
	rMin = 2**(1./6.)
	m = 1*(math.pi*pow(rMin,3))/(15*T)
	m += (4./3)*math.pi*(1/T)*((1/(7*pow(rc,7.))-1/(1*pow(rc,1)))-(1/(7*pow(rMin,7.))-1/(1*pow(rMin,1))))
	return m

def fmu_id(rho,M):
	rhoM = rho/M
	return math.log(rhoM)

def fp_id(rho,M):
	rhoM = rho/M
	return rhoM

def fmu_hs(rho):
	eta = rho*math.pi/6
	return eta*(8-9*eta+3*eta*eta)/(pow(1-eta,3.))

def fp_hs(rho):
	eta = rho*math.pi/6
	return rho*(4*eta-2*eta*eta)/(pow(1-eta,3.))

def fmu_disp(rho,a):
	return -2*a*rho

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

def get_p(rho,a,Ma,K,T_assoc,M):
	p = fp_id(rho,M) + (fp_hs(rho) + fp_disp(rho,a) + fp_assoc(rho,Ma,K,T_assoc) + fp_chain(rho,M))
	return p


def find_peaks_mu(T,a,Ma,K,T_assoc,M):
	denVec = np.arange(0.0001,1.5,0.0001)

	test_solVec = []
	for i in range (1,len(denVec)-1):
		f = get_mu(denVec[i],a,Ma,K,T_assoc,M)
		f_before = get_mu(denVec[i-1],a,Ma,K,T_assoc,M)
		f_after = get_mu(denVec[i+1],a,Ma,K,T_assoc,M)

		if (f-f_before)*(f_after-f)<0:
			test_solVec.append(denVec[i])

	mu_max = get_mu(test_solVec[0],a,Ma,K,T_assoc,M)
	mu0 = get_mu(denVec[0],a,Ma,K,T_assoc,M)
	mu_min = get_mu(test_solVec[1],a,Ma,K,T_assoc,M)
	if mu_min<mu0:
		mu_min = mu0

	return [mu_max,mu_min]


def find_roots_mu(mu_test,T,a,Ma,K,T_assoc,M):
	tol_ref = 1e-6
	denVec = np.arange(0.0001,2.0,0.0001)

	test_solVec = []
	for i in range (1,len(denVec)):
		f = get_mu(denVec[i],a,Ma,K,T_assoc,M) - mu_test
		f_before = get_mu(denVec[i-1],a,Ma,K,T_assoc,M) - mu_test

		if f>0 and f_before<0:
			test_solVec.append(denVec[i-1])
			test_solVec.append(denVec[i])
	
	tol = 1.
	step = 0
	while (tol>tol_ref):
		if (step == 0):
			sol1 = test_solVec[0]
			sol2 = test_solVec[1]

		sol3 = (sol1+sol2)/2

		f1 = get_mu(sol1,a,Ma,K,T_assoc,M) - mu_test
		f2 = get_mu(sol2,a,Ma,K,T_assoc,M) - mu_test
		f3 = get_mu(sol3,a,Ma,K,T_assoc,M) - mu_test


		if f3>0 and f1<0:
			sol2 = sol3
		elif f2>0 and f3<0:
			sol1 = sol3
		else:
			print ('Something is wrong..')
			exit(1)
		tol = abs(sol2-sol1)
		step += 1
	rho_v = (sol1+sol2)*0.5

	tol = 1.
	step = 0
	while (tol>tol_ref):
		if (step == 0):
			sol1 = test_solVec[2]
			sol2 = test_solVec[3]

		sol3 = (sol1+sol2)/2

		f1 = get_mu(sol1,a,Ma,K,T_assoc,M) - mu_test
		f2 = get_mu(sol2,a,Ma,K,T_assoc,M) - mu_test
		f3 = get_mu(sol3,a,Ma,K,T_assoc,M) - mu_test



		if f3>0 and f1<0:
			sol2 = sol3
		elif f2>0 and f3<0:
			sol1 = sol3
		else:
			print ('Something is wrong..')
			exit(1)
		tol = abs(sol2-sol1)
		step += 1
	rho_l = (sol1+sol2)*0.5

	return [rho_v, rho_l]



Tc = 1.626415
#T_vec = np.arange(Tc*0.7,Tc+0.001,0.1)
T_vec = np.array([0.7*Tc,0.75*Tc,0.8*Tc,0.85*Tc,0.9*Tc,0.95*Tc,0.995*Tc])
T_vec = np.array([0.7*Tc])

M = 1.
Ma = 4.
K = 1.4849*1e-4
T_assoc = 1/5.

rc = 5.0
rMin = 2**(1./6.)
dr = 0.01

o1filename = 'coex_eps'+str(1/T_assoc)+'_M'+str(M)+'_Ma'+str(Ma)+'.xvg'
o2filename = 'mu0_eps'+str(1/T_assoc)+'_M'+str(M)+'_Ma'+str(Ma)+'.xvg'
o3filename = 'm_eps'+str(1/T_assoc)+'_M'+str(M)+'_Ma'+str(Ma)+'.xvg'
o4filename = 'a_eps'+str(1/T_assoc)+'_M'+str(M)+'_Ma'+str(Ma)+'.xvg'



oMat = []
o2Mat = []
o3Mat = []
o4Mat = []
for i in range (len(T_vec)):
	T = T_vec[i]
	a = get_a(T,rc)
	m = get_m(T,rc)

	# normalizeation T_assoc relative to the T
	T_assoc_norm = T_assoc


	print "T = "+str(T)

	# find mu_coex and p_coex
	mu_max, mu_min = find_peaks_mu(T,a,Ma,K,T_assoc_norm,M)
	print "mu_max, mu_min = "+str(mu_max)+" "+str(mu_min)


	#mu_test = mu_min + random.random()*(mu_max-mu_min)
	mu_test = (mu_max+mu_min)*0.5
	print "mu_test = "+str(mu_test)

	rho_v, rho_l = find_roots_mu(mu_test,T,a,Ma,K,T_assoc_norm,M)

	pv = get_p(rho_v,a,Ma,K,T_assoc_norm,M)
	pl = get_p(rho_l,a,Ma,K,T_assoc_norm,M)
	#print pv, pl

	if pv > pl:
		mu1 = mu_test
		mu2 = mu_max
	elif pv < pl:
		mu1 = mu_min
		mu2 = mu_test

	tol_ref = 1e-8
	tol = 1
	pl_before = 100.
	while (tol>tol_ref):
		mu3 = (mu1+mu2)*0.5
		rho_v, rho_l = find_roots_mu(mu3,T,a,Ma,K,T_assoc_norm,M)

		pv = get_p(rho_v,a,Ma,K,T_assoc_norm,M)
		pl = get_p(rho_l,a,Ma,K,T_assoc_norm,M)

		if pv > pl:
			mu1 = mu3
		elif pv < pl:
			mu2 = mu3
		tol = abs(mu2-mu1)

		print rho_v, rho_l, pv, pl, pl_before


	mu_coex = (mu1+mu2)/2
	p_coex = (pl+pv)/2

	print rho_v, rho_l, mu_coex, p_coex
	print fmu_id(rho_v,M) + M*(fmu_hs(rho_v) + fmu_disp(rho_v,a) + fmu_chain(rho_v,M) + fmu_assoc(rho_v,Ma,K,T_assoc_norm)), fmu_id(rho_l,M) + M*(fmu_hs(rho_l) + fmu_disp(rho_l,a) + fmu_chain(rho_l,M) + fmu_assoc(rho_l,Ma,K,T_assoc_norm))
	print fp_id(rho_v,M) + (fp_hs(rho_v) + fp_disp(rho_v,a) + fp_chain(rho_v,M) + fp_assoc(rho_v,Ma,K,T_assoc_norm)) ,fp_id(rho_l,M) + (fp_hs(rho_l) + fp_disp(rho_l,a) + fp_chain(rho_l,M) + fp_assoc(rho_l,Ma,K,T_assoc_norm))



# we find mu_coex, rho_v, rho_l
sVec = [0.01,0.05,0.1,0.15,0.175,0.2]

deltamu = 1e-5
for i in range (len(sVec)):
	try:
		s = sVec[i]
		mu = mu_coex*(1-s)
		#important
		mu = mu - deltamu
		rho_v, rho_l = find_roots_mu(mu,T,a,Ma,K,T_assoc,M)

		pv = get_p(rho_v,a,Ma,K,T_assoc,M)
		pl = get_p(rho_l,a,Ma,K,T_assoc,M)

		print s, fp_id(rho_l,M), fp_hs(rho_l), fp_disp(rho_l,a), pl, pv


		oMat.append([T, s, mu_coex, mu, rho_l, rho_v, pl, pv])
	except:
		pass

o1filename = 'uncoex_eps'+str(1/T_assoc)+'_M'+str(M)+'_Ma'+str(Ma)+'_T'+str(T)+'_dmu'+str(deltamu)+'.xvg'

np.savetxt(o1filename,oMat,fmt='%5.10f')


'''
	if abs(pl - pl_before) < 1e-8:
		print "change the iteration direction"
		if pv < pl:
			mu1 = mu_test
			mu2 = mu_max
		elif pv > pl:
			mu1 = mu_min
			mu2 = mu_test

		while (tol>tol_ref):
			mu3 = (mu1+mu2)*0.5
			rho_v, rho_l = find_roots_mu(T,a,mu3)
			pv = fp(rho_v,T,a)
			pl = fp(rho_l,T,a)
			if pv < pl:
				mu1 = mu3
			elif pv > pl:
				mu2 = mu3
			tol = abs(pv-pl)

			print rho_v, rho_l, pv, pl
'''

