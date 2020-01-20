#!/usr/bin/env python 
from matplotlib import pyplot as plt
import matplotlib.mpl as mpl
import matplotlib.cm as cm
import numpy as np
import numpy.linalg as la
import cPickle as pickle
import ndfit
import math
import sys

## For timing siulations 
import timeit

## For Plotting
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':20})
rc('text', usetex=True)
rc('xtick', labelsize=20) 
rc('ytick', labelsize=20)
rc('legend', fontsize=20) 
#rc({'legend.numpoints':1 ,'legend.scatterpoints':1})

## Set max recursion depth
depth = int(1000)
sys.setrecursionlimit(depth)

## A class with two methods which to calclates the DFT
## and IDFT matricies. One passes to fft the frequency
## one is interested in, and the number of harmoics to
## generate

# WARNING: This script contains known errors (archive only)

class fft:

	def __init__(self, freq, n):

		## Basic stuff
		self.n = n
		omega = 2*math.pi*freq
		tau = 1./freq
		
		## Get the vectors of frequencies and times
		self.omega = [i*omega for i in range(n)]
		self.tau = [tau*i/n for i in range(n)]

		## Calculate and save the dft and idft
		self.dft  = self.dft()
		self.idft = self.idft()

	## Note that the DFT and IDFT transformations the inverse is the adjoint matrix
	## up to the constant n. (The number of Harmonics) 
	def dft(self):
		dft = np.zeros(shape=(self.n,self.n), dtype='complex')
		for i in range(self.n):
			for j in range(self.n):
				dft[i][j] = np.exp(np.complex(0,-self.omega[i]*self.tau[j]))
		return dft


	def idft(self):
		idft = np.zeros(shape=(self.n,self.n), dtype='complex')
		for i in range(self.n):
			for j in range(self.n):
				idft[i][j] = np.exp(np.complex(0,self.omega[j]*self.tau[i]))
		return idft


	def fft(self,v): return np.dot(self.dft,v)
	def ifft(self,v): return np.dot(self.idft,v)
		

## A class which forms the dual vector of a vector f
## Using the DFT/IDFT method. Mode specifys whether
## The vector passed is in the time or frequency
## domain. Every dual will also have a vector in
## the time domain associated with it. 
class dual:

	def __init__(self, v, mode,f=None):

		if f == None:
			self._fft = fft(1./len(v),len(v))
			self.tau  = self._fft.tau 
		else:
			self._fft = fft(f,len(v))
			self.tau  = self._fft.tau 
		if mode == 't':
			self.t = v
			self.f = np.dot(self._fft.dft,v)
		if mode == 'f':
			self.t = np.dot(self._fft.idft,v)
			self.f = v


## Make a diode class. This stores our diode model. This is the Rizoli model.
## Its effect is to weaken the non-linearity in the diode		
class diode:

	## Allow the resistnace to be changed
	def __init__(self):

		## Constants of the problem
		self.vstop = complex(0.85,0)
		self.vt,self.n,self.Is,self.Cj,self.phi,self.g = 0.026,1.3,1e-11,0.15e-12,1.3,0.5

	## The diode capacitance
	def diodec(self,vd):
		
		return self.Cj*np.power(1-(vd/self.phi), -self.g)

	## The diode current
	def diodei(self,vd):
		
		if vd<self.vstop:
		
			return (self.Is)*(np.exp(vd/self.n/self.vt)) - self.Is
		
		else:
		
			return (self.Is)*(np.exp(self.vstop/self.n/self.vt))*(1 + ((vd-self.vstop)/self.n/self.vt))-self.Is

	## The condutance
	def diodeg(self,vd):
		
		if vd<self.vstop:
		
			return (self.Is/self.n/self.vt)*(np.exp(vd/self.n/self.vt))
		
		else:
			
			return (self.Is)*(np.exp(self.vstop/self.n/self.vt))*(vd/self.n/self.vt)

## This is diode class which is identical to the above with the exception
## that the polarity of the diode is reversed. This is neede for cross diode
## simulations.
class diodex:

	## Allow the resistnace to be changed
	def __init__(self):

		## Constants of the problem
		self.vstop = complex(-0.85,0)
		self.vt,self.n,self.Is,self.Cj,self.phi,self.g = 0.026,1.3,1e-11,0.15e-12,1.3,0.5

	## The diode capacitance
	def diodec(self,vd):
		
		return self.Cj*np.power(1-(-vd/self.phi),-self.g)
		
	## The diode current
	def diodei(self,vd):
		if vd>self.vstop:
			
			return -(self.Is*np.exp(-1*vd/self.n/self.vt)-self.Is)
		else:
			
			return -(self.Is*np.exp(-1*self.vstop/self.n/self.vt)*(1-(vd-self.vstop)/self.n/self.vt) - self.Is)

	## The condutance
	def diodeg(self,vd):

		if vd>self.vstop:
			
			return (self.Is/self.n/self.vt)*(np.exp(-vd/self.n/self.vt))

		else:
			
			return self.Is*(np.exp(-self.vstop/self.n/self.vt))*(-vd/self.n/self.vt)

def harmonicBalance(_V, _Vs, freq, npoints, CONV=5.0e-4, counter=0, J=None):

	# Set convergence and recursion counter
	_CONV = CONV
	counter+=1

	## Establish the number of harmonics and initialize
	## the diode model
	d,x = diode(),diodex()
	omega = 2*math.pi*freq
	
	G,C1,C2,L1,R= 0.02,0.15e-12,0.15e-12,1.64e-9,1.5
	## This is the admittance from port 1 to port 1. This is
	## our shunt element 
	Y11 = np.zeros(shape=(npoints,npoints), dtype='complex')
	for i in range(npoints):
		for j in range(npoints):
			if i == j: 
				Y11[i][j] = complex(2*G,j*omega*C1/(1-(j*j*omega*omega*L1*C1)+complex(0,j*omega*R*C1))+i*omega*C2) 
			else:
				Y11[i][j] = 0.0

	## Calculate the linear admittance matrix. This is Y12.
	## This is our series element (SOURCE ADMITTANCE 50OHM)
	Y12 = np.zeros(shape=(npoints,npoints), dtype='complex')
	for i in range(npoints):
		for j in range(npoints):
			if i == j: 
				Y12[i][j] = complex(-G,0)
			else:
				Y12[i][j] = 0.0

	## Diagonal matrix of harmonics
	OO = np.zeros(shape=(npoints,npoints), dtype='complex')
	for i in range(npoints):
		for j in range(npoints):
			if i == j: 
				OO[i][j] = complex(0.0,i*omega)
			else:
				OO[i][j] = 0.0

	## Calculate the operator F for the non-linearity
	## and its dual (time domain). This is actually
	## a current

	## Forward Diode implementation 
	Fv = [d.diodei(_V.t[i]) for i in range(npoints)]

	## Reverse Diode implementation
	#Fv = [x.diodei(_V.t[i]) for i in range(npoints)]

	## Cross Diode implementation
	#Fv = [(x.diodei(_V.t[i])+d.diodei(_V.t[i])) for i in range(npoints)]
	_VFv = dual(Fv,'t',freq)

	## Calculate the error function and its dual
	## (frequency domain)
	ef  = np.dot(Y11,_V.f) + np.dot(Y12,_Vs.f) + _VFv.f
	_Ef = dual(ef,'f',freq)

	#############################
	#	 RECURSION BLOCK	 #
	#############################
	print max(np.abs(ef))/npoints
	if max(np.abs(ef))/npoints < _CONV:
		return _V,J

	elif counter > depth-20:
		print "Hit recursion depth. Returning something"
		return _V,J
	
	else:
		## Calculate the Jacobian. Recall that it is
		## diagonal because we have a LTI system.
		## (time domain). Set epsilon as a throttle on
		## stepping the delta vector.
		epsilon = 1.1
		
		## Resistive term in Jacobian
		DG = np.zeros(shape=(npoints,npoints), dtype='complex')
		for i in range(npoints):
			for j in range(npoints):
				if i == j: 
					DG[i][j] = d.diodeg(_V.t[i]) + x.diodeg(_V.t[i])
				else:
					DG[i][j] = 0.0

		## Capacitive term in Jacobian
		DC = np.zeros(shape=(npoints,npoints), dtype='complex')
		for i in range(npoints):
			for j in range(npoints):
				if i == j: 
					DC[i][j] = d.diodec(_V.t[i]) + x.diodec(_V.t[i])
				else:
					DC[i][j] = 0.0

		## Get our fft from the fft class and calculate change
		## in the error function (frequency domain)			
		T = fft(freq,npoints)
		J = Y11 + np.dot(T.dft,np.dot(DG,T.idft)) + np.dot(OO,np.dot(T.dft,np.dot(DC,T.idft)))
	 
		## Invert this matrix and multiply it by the error
		## to obtain change in voltage (frequency domain).
		## add it to original vector V and calculate its dual.
		## (frequency domain).
		DV  = np.dot(la.inv(J),ef)
		_V1 = dual(_V.f - epsilon*DV,'f',freq)

		## iterate 
		return harmonicBalance(_V1, _Vs, freq, npoints,_CONV, counter, J)

if __name__ == "__main__":

	###############################
	# LARGE SIGNAL HB		
	#
	if 1 == 1:
		## Establish the number of harmonics and initialize the diode model
		npoints, freq, VOLTAGE, CONV = 50, 5e9, [0.75,1.5,3.0], [1.0e-5, 1.94e-4, 4.95e-4]

		## Set initialization vector. 
		v = np.zeros([npoints], dtype="complex")
		v[0], v[1], v[3] = 0.0, 0.0, 0.0
		_V = dual(v,'f',freq)

		# For loop for three different drive voltages
		RESULT,SOURCE = [],[]
		IT, GT, CT = [],[],[]
		for i,_VOLTAGE in enumerate(VOLTAGE):
		
			## Set a source vector in the frequency domain
			vs = np.zeros([npoints], dtype="complex")
			vs[0], vs[1], vs[3] = 0, _VOLTAGE*complex(0,-1), 0.0
			_Vs = dual(vs,'f',freq) 

			## Run the iterations and return the dual which is a solution 
			print "-----------------------------------------------"
			V,Jacobian = harmonicBalance(_V, _Vs, freq, npoints, CONV[i])
			RESULT.append(V), SOURCE.append(_Vs)
		
			## Get the other waveforms
			f,r= diode(),diodex()

			## For Regular Diode
			#IT.append([f.diodei(V.t[i]) for i in range(npoints)])
			#GT.append([f.diodeg(V.t[i]) for i in range(npoints)])
			#CT.append([f.diodec(V.t[i]) for i in range(npoints)])

			## For Cross Diode
			IT.append([f.diodei(V.t[i])+r.diodei(V.t[i]) for i in range(npoints)])
			GT.append([f.diodeg(V.t[i])+r.diodeg(V.t[i]) for i in range(npoints)])
			CT.append([f.diodec(V.t[i])+r.diodec(V.t[i]) for i in range(npoints)])

		#DATA  = {"RESULT": RESULT, "CT":CT, "GT":GT}
		#pickle.dump(DATA, open('DATA.pkl', 'wb'))

	################################
	# FOR DC AND RF SIMULATIONS
	#
	if 1 == 1:

		fig = plt.figure(1)
		ax1 = plt.subplot(211)
		ax1.set_title("Voltage Waveform")
		for R in RESULT:
			ax1.plot(R.t,'-')
		ax1.set_ylim(-3,3)
		ax1.set_ylabel("Voltage $(V)$")
		ax1.set_xlim(0,44)
		ax2 = plt.subplot(212)
		for S in SOURCE:
			ax2.plot(S.t,'--')
		ax2.set_xlabel("Time Step $(n)$")
		ax2.set_ylabel("Voltage $(V)$")
		ax1.set_ylim(-1,1)
		ax2.set_xlim(0,44)

		fig = plt.figure(2)
		ax1 = plt.subplot(111)
		ax1.set_title("Conversion Loss")
		hlist = []
		for i,R in enumerate(RESULT):
			hlist.append(ax1.plot(-20*np.log10(VOLTAGE[i]/np.abs(R.f)),'o'))
		ax1.set_xlim(0,9)
		ax1.set_ylim(-45,-5)
		ax1.legend(hlist,['$0.75V$','$1.50V$','$3.00V$'],numpoints=1) 
		ax1.set_xlabel("Harmonic Number $(n)$")
		ax1.set_ylabel("Conversion Loss $-10\log(P_{RF}/P_{IF})$")
			
		fig = plt.figure(3)
		ax1 = plt.subplot(111)
		for It in IT: 
			ax1.semilogy(np.abs(It))
		ax1.set_xlabel("Time Step $(n)$")
		ax1.set_ylabel("Current Waveform $|I|$")
		ax2 = ax1.twinx()
		for Gt in GT:
			ax2.semilogy(np.abs(Gt),'--')
		ax2.set_ylabel("Conductance Waveform $G$")
		ax1.set_xlim(0,44)

		fig = plt.figure(4)
		ax1 = plt.subplot(111)
		ax1.set_xlabel("Time Step $(n)$")
		ax1.set_ylabel("Capacitance Waveform $C$")
		for Ct in CT:
			plt.plot(Ct)
		ax1.set_xlim(0,44)

		plt.show()
