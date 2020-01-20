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

#####################################
#  Methods for Conversion Matricies #
#####################################
def CONVERSION(num, fft):
	# Construct the conversion matrix
	size = (2*num+1)
	CMATRIX = np.zeros((size,size), dtype="complex")
	for i,v in enumerate(fft[0:size]):
		if i == 0: CMATRIX+=np.diag([v]*(size))	   
		else: CMATRIX+= (np.diag([np.conj(v)]*(size-i),i) + np.diag([v]*(size-i),-i))	 
	#print np.round(CMATRIX,3)
	return CMATRIX

def OMEGA(num, freq, sub):
	# Construct the omega matrix
	size = (2*num+1)
	V_OMEGA = [i*freq + sub for i in range(-num, num+1)] 
	M_OMEGA = np.zeros((size,size),dtype="complex")
	for i,f in enumerate(V_OMEGA): M_OMEGA[i][i] = complex(0,f)
	return M_OMEGA,V_OMEGA
	
def PASSIVES(vomega):
	G,C1,C2,L1,R = 0.02,0.15e-12,0.15e-12,1.64e-9,1.5
	Y = np.zeros((len(vomega),len(vomega)),dtype="complex")
	for i,omega in enumerate(vomega):
		Y[i][i] = np.complex(0,omega*C1/(1-(omega*omega*L1*C1)+complex(0,omega*R*C1)) + omega*C2) 
	return Y

def DIAGONAL(vomega,zl=50):
	D = np.zeros((len(vomega),len(vomega)),dtype="complex")
	for i,omega in enumerate(vomega): D[i][i] = complex(zl,0) 
	return D
	
def BUILDTWOPORT(y11, y12, y21, y22):

	if y11 is not None: TWOPORT = np.zeros((2*len(y11),2*len(y11)))
	if y22 is not None: TWOPORT = np.zeros((2*len(y22),2*len(y22)))
		
	if y11 is not None:
		for i,value in np.ndenumerate(y11): TWOPORT[i[0],i[1]] = value
	if y12 is not None:
		for i,value in np.ndenumerate(y12): TWOPORT[i[0],i[1]+len(y22)] = value
	if y21 is not None:
		for i,value in np.ndenumerate(y21): TWOPORT[i[0]+len(y22),i[1]] = value
	if y22 is not None:
		for i,value in np.ndenumerate(y22): TWOPORT[i[0]+len(y22),i[1]+len(y22)] = value
	return TWOPORT

if __name__ == "__main__":

	####################################################
	##   CONVERSION MATRICIES FOR TWO TONE ANALYSIS   ##
	####################################################
	
	# Select the number of harmonics to simulate to.
	# We would like data up to 16GHz (at least)
	num = 3

	# Obtain the data from the single diode simulation
	DATA, freq = pickle.load(open('DATA.pkl', 'rb')), 5e9
	RESULT,CT, GT = DATA["RESULT"],DATA["CT"],DATA["GT"]
		
	# Calculate the fft of conductance and capacitance waveforms
	_GT = [np.fft.fft(g)/len(g) for g in GT]
	_CT = [np.fft.fft(c)/len(c) for c in CT]

	# Set up the omega vector and frequency
	OMEGA_M, OMEGA_V = OMEGA(num , 5e9, 1e9) 

	hlist = []
	for i in range(len(_GT)):

		# Calculate conversion matricies for the diode and passives
		G_CONVERSION = CONVERSION(num ,_GT[i])
		C_CONVERSION = np.dot(OMEGA_M, CONVERSION(num,_CT[i]))
		D_CONVERSION = G_CONVERSION + C_CONVERSION 
		Y_CONVERSION = PASSIVES(OMEGA_V)

		# Calculate the conversion matrix of the total network
		# (e.g. a bunch of shunt elements)
		NETWORK = D_CONVERSION + Y_CONVERSION 
	
		# Calculate SOURCE and LOAD conversion matricies (50 Ohm)
		_SOURCE, _LOAD, _ONES = DIAGONAL(OMEGA_V), DIAGONAL(OMEGA_V), DIAGONAL(OMEGA_V,1e9)  

		# Build the source and load components of the
		# twoport conversion matrix
		Z11,Z12,Z21,Z22 = _SOURCE, None, None, _LOAD
		ZS  =  BUILDTWOPORT(Z11,Z12,Z21,Z22)
		
		# Build the network and twoport conversion matrix
		Y11,Y12,Y21,Y22 = _ONES, -1*_ONES, -1*_ONES, NETWORK + _ONES
		YNET =  BUILDTWOPORT(Y11,Y12,Y21,Y22)
		YF   =  la.inv(ZS + la.inv(YNET))

		# Select the frequencey one would like to extract
		omega_r, size = 1e9, len(Z11) 

		GAIN = []
		for omega_s in OMEGA_V :
			IR, IS = OMEGA_V.index(omega_r), OMEGA_V.index(omega_s) 
	
			# Get the y-parameters of the value in question
			y11 = YF[IR][IR]
			y12 = YF[IR][IS+size]
			y21 = YF[IS+size][IR]
			y22 = YF[IS+size][IS+size]

			# Calculate the conversion loss
			YS,YL = 0.02,0.02
			g = 4*YS*YL*(np.abs(y21)**2)/( np.abs((y11 + YS)*(y22 + YL) - y12*y21)**2)
			GAIN.append(10*np.log10(g))

		hlist.append(plt.plot([w/1e9 for w in OMEGA_V], GAIN, 'o'))

	plt.legend(hlist,["$0.75V$","$1.50V$","$3.00V$"],loc=1,numpoints=1)
	plt.title("Up Conversion")
	plt.xlabel("Frequency $(GHz)$")
	plt.ylabel("Conversion Loss $(dB)$")
	plt.ylim(-35,0)
	plt.show()
