import numpy as np
import sys
import os

Ri = float(sys.argv[1])
R = 40.

mu1 = -2.44
mu2 = -2.491

deriv1 = 0
deriv2 = 0

while deriv1*deriv2 >= 0:
	o1FileName = 'denout_Ri'+str(Ri)+'_R'+str(R)+'_mu'+str(mu1)+'.xvg'
	o2FileName = 'denout_Ri'+str(Ri)+'_R'+str(R)+'_mu'+str(mu2)+'.xvg'
	om1FileName = 'monitor_Ri'+str(Ri)+'_R'+str(R)+'_mu'+str(mu1)+'.xvg'
	om2FileName = 'monitor_Ri'+str(Ri)+'_R'+str(R)+'_mu'+str(mu2)+'.xvg'

	os.system('./qdft input_test '+o1FileName+' '+om1FileName+' '+str(mu1)+' '+str(R)+' '+str(Ri)+' 0')
	os.system('./qdft input_test '+o2FileName+' '+om2FileName+' '+str(mu2)+' '+str(R)+' '+str(Ri)+' 0')

	o1Mat = np.loadtxt(om1FileName)
	o2Mat = np.loadtxt(om2FileName)
	o1Vec = np.transpose(o1Mat)[1]
	o2Vec = np.transpose(o2Mat)[1]

	for i in range (len(o1Vec)):
		deriv1 = o1Vec[len(o1Vec)-1]-o1Vec[len(o1Vec)-2-i]
		if abs(deriv1) > 2.5:
			break
	for i in range (len(o2Vec)):
		deriv2 = o2Vec[len(o2Vec)-1]-o2Vec[len(o2Vec)-2-i]
		if abs(deriv2) > 2.5:
			break
	print deriv1, deriv2

	while deriv2 >= 0:
		mu3 = mu2*1.02
		o3FileName = 'denout_Ri'+str(Ri)+'_R'+str(R)+'_mu'+str(mu3)+'.xvg'
		om3FileName = 'monitor_Ri'+str(Ri)+'_R'+str(R)+'_mu'+str(mu3)+'.xvg'
		os.system('./qdft input_test '+o3FileName+' '+om3FileName+' '+str(mu3)+' '+str(R)+' '+str(Ri)+' 40')
		o3Mat = np.loadtxt(om3FileName)
		o3Vec = np.transpose(o3Mat)[2]
		for i in range (len(o3Vec)):
			deriv3 = o3Vec[len(o3Vec)-1]-o3Vec[len(o3Vec)-2-i]
			if abs(deriv3) > 3:
				break

		if deriv1*deriv3 < 0:
			mu2 = mu3
		elif deriv2*deriv3 < 0:
			mu1 = mu3

		os.system('rm -f '+o3FileName)
		os.system('rm -f '+om3FileName)

	while deriv1 < 0:
		mu3 = mu1*0.98
		o3FileName = 'denout_mu'+str(Ri)+'_R'+str(R)+'_mu'+str(mu3)+'.xvg'
		om3FileName = 'monitor_mu'+str(Ri)+'_R'+str(R)+'_mu'+str(mu3)+'.xvg'
		os.system('./qdft input_test '+o3FileName+' '+om3FileName+' '+str(mu3)+' '+str(R)+' '+str(Ri)+' 40')
		o3Mat = np.loadtxt(om3FileName)
		o3Vec = np.transpose(o3Mat)[2]
		for i in range (len(o3Vec)):
			deriv3 = o3Vec[len(o3Vec)-1]-o3Vec[len(o3Vec)-2-i]
			if abs(deriv3) > 3:
				break

		if deriv1*deriv3 < 0:
			mu2 = mu3
		elif deriv2*deriv3 < 0:
			mu1 = mu3

		os.system('rm -f '+o3FileName)
		os.system('rm -f '+om3FileName)


	tol = abs(mu1-mu2)

	os.system('rm -f '+o1FileName)
	os.system('rm -f '+o2FileName)
	os.system('rm -f '+om1FileName)
	os.system('rm -f '+om2FileName)
os.system('rm -f *.energy')


while tol > 1E-7:
	mu3 = (mu1+mu2)/2
	o3FileName = 'denout_mu'+str(Ri)+'_R'+str(R)+'_mu'+str(mu3)+'.xvg'
	om3FileName = 'monitor_mu'+str(Ri)+'_R'+str(R)+'_mu'+str(mu3)+'.xvg'
	os.system('./qdft input_test '+o3FileName+' '+om3FileName+' '+str(mu3)+' '+str(R)+' '+str(Ri)+' 0')
	o3Mat = np.loadtxt(om3FileName)
	o3Vec = np.transpose(o3Mat)[1]
	for i in range (len(o3Vec)):
		deriv3 = o3Vec[len(o3Vec)-1]-o3Vec[len(o3Vec)-2-i]
		if abs(deriv3) > 3:
			break

	if deriv1*deriv3 < 0:
		mu2 = mu3
	elif deriv2*deriv3 < 0:
		mu1 = mu3

	tol = abs(mu1-mu2)

	print "mu1 = "+str(mu1)+", mu2 = "+str(mu2)+", tol = "+str(tol)

	os.system('rm -f '+o3FileName)
	os.system('rm -f '+om3FileName)
os.system('rm -f *.energy')

print "mu* = "+str((mu1+mu2)/2)

mu4 = (mu1+mu2)/2
o4FileName = 'denout_Ri'+str(Ri)+'_R'+str(R)+'_mu'+str(mu4)+'.xvg'
om4FileName = 'monitor_Ri'+str(Ri)+'_R'+str(R)+'_mu'+str(mu4)+'.xvg'
os.system('./qdft input_run '+o4FileName+' '+om4FileName+' '+str(mu4)+' '+str(R)+' '+str(Ri)+' 0')



