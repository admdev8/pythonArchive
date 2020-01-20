#!/usr/bin/env python
import math
import numpy as np
import random as r
import matplotlib.pyplot as plt
import multiprocessing as mp
from collections import defaultdict

#################################################################
#-------------- GaAs Monte Carlo Simulation Code ---------------#
#################################################################
# Intravalley emission and absorption of optical and acoustic #
# phonons and intervalley emission and absorption of phonons. #
#################################################################
#############################################
# A container class for universal constants #
#############################################
class Universe(object):

	def __init__(self):

		# UNIVERSAL CONSTANTS
		self.e = 1.602e-19
		#
		self.pi = 3.1415927
		#
		self.kT = 0.026
		#
		self.hbar = 6.5821e-16
		#
		self.ep = 8.8548e-14*self.e
		#
		self.c = 3e10
		#
		self.me = 0.511e6/(self.c**2)
		#
		#Fundamental Charge (C)
		#pi (unitless)
		#Thermal Energy (eV)
		#Planck Constant (eV/s)
		#permitivvity (eV/V^2)
		#Speed of Light (cm/s)
		#Electron Mass (eV/c^2)



#####################################
# A container class for GaAs values #
#####################################
class GaAs(Universe):

	def __init__(self):

		Universe.__init__(self);

		# STATIC CONSTATS
		self.rho = 3.0123e33/(self.c**2)
		self.nu = 5.22e5
		self.epI = 10.82
		self.epN= 12.53

		# PHONONS
		self.wOP= 5.37e13
		self.wE= 4.54e13

		# DEFORMATION POTENTIALS
		self.Da= 7.0
		self.De= 1e9

		# VALLEY EFFECTIVE MASSES
		self.mc= 0.067*self.me
		self.ms= 0.350*self.me

		# VALLEY SEPERATION
		self.D= 0.36
		self.Z= 3.0
		
		#Density (eV/c^2cm^3)
		#Speed of Sound (cm/s)
		#High Frequency Dielectric Const
		#Static Dielectric Const
		# Polar OP Phonon (rad/s)
		# Eq Intervalley Phonon (rad/s)
		# Acoustic (eV)
		# Eq Intervalley (eV)
		# Central Valley (g)
		# Satellite Valley (g)
		# Central to Sattelite (eV)
		# Number of equivalent valleys (#)


#########################################
# VALLEY superclass
#
#########################################
# This class contains methods which calculate all of the
# scattering lifetimes. These methods will be called from
# two valley subclasses (one for G and one for L).
class Valley(GaAs):

	def __init__(self, V):

		GaAs.__init__(self);

		if V=="G": 
			self.M = self.mc
		
		if V=="L": 
			self.M = self.ms

	def AcousticPhonons(self, E):
		
		a = np.power(2*(self.M),1.5)*self.kT*(self.Da**2)
		b = 4*self.pi*self.rho*(self.nu**2)*(self.hbar**4)
		
		result = (a/b)*np.sqrt(E)
		
		for i,r in enumerate(result):
			if np.isnan(r): result[i] = 0.0
		return result

	def OpticalPhonons(self, E, mode):
		Nop = 1.0/(np.exp(self.wOP*self.hbar/self.kT) - 1)
		
		if mode=="A": 
			EO = E+(self.hbar*self.wOP) 
			Nop+=0
		
		if mode=="E": 
			EO = E-(self.hbar*self.wOP) 
			Nop+=1
		
		a = ((Nop*(self.e**2)*np.sqrt(self.M)*self.wOP))
		b = (np.sqrt(2)*self.hbar*(4*self.pi*self.ep))
		c = ((1.0/self.epI) - (1.0/self.epN))
		factor = (a*c)/b
	
		arg = abs((np.sqrt(E) + np.sqrt(EO))/(np.sqrt(E) - np.sqrt(EO)))
		result = factor*np.log(arg)/np.sqrt(E)
	
		for i,r in enumerate(result):
			if np.isnan(r): result[i] = 0.0

		return result

	def CentralToSatellite(self, E, mode):
	
		Nop = 1.0/(np.exp(self.wOP*self.hbar/self.kT) - 1)
		
		if mode=="A": 
			EO = E+(self.hbar*self.wOP)-self.D
			Nop+=0
		
		if mode=="E": 
			EO = E-(self.hbar*self.wOP)-self.D 
			Nop+=1
		
		a = np.power(3*self.ms,1.5)*(self.De**2)*Nop
		b = np.sqrt(2)*self.pi*self.rho*self.wE*(self.hbar**3)
		result = (a/b)*np.sqrt(EO)
		
		for i,r in enumerate(result):
			if np.isnan(r): result[i] = 0.0
		return result

	def SatelliteToCentral(self, E, mode):

		Nop = 1.0/(np.exp(self.wOP*self.hbar/self.kT) - 1)
	
		if mode=="A": 
			EO = E+(self.hbar*self.wOP)+self.D 
			Nop+=0
	

		if mode=="E": 
			EO = E-(self.hbar*self.wOP)+self.D
			Nop+=1
	
		a = np.power(self.mc,1.5)*(self.De**2)*Nop
		b = np.sqrt(2)*self.pi*self.rho*self.wE*(self.hbar**3)
		result = (a/b)*np.sqrt(EO)
	
		for i,r in enumerate(result):
			if np.isnan(r): result[i] = 0.0
	
		return result
	
	def SatelliteToSatellite(self,E,mode):
	
		Nop = 1.0/(np.exp(self.wOP*self.hbar/self.kT) - 1)
		if mode=="A": 
			EO = E+(self.hbar*self.wOP) 
			Nop+=0

		if mode=="E": 
			EO = E-(self.hbar*self.wOP)
			Nop+=1
		
		a = (self.Z - 1)*np.power(self.ms,1.5)*(self.De**2)*Nop
		b = np.sqrt(2)*self.pi*self.rho*self.wE*(self.hbar**3)
		result = (a/b)*np.sqrt(E)
		
		for i,r in enumerate(result):
			if np.isnan(r): result[i] = 0.0
		
		return result


################################################
#VALLEY subclasses
#
################################################
class g_Valley(Valley):

	def __init__(self, E):
	
		super(g_Valley,self).__init__("G")
		self.E = E
		self.AC = super(g_Valley, self).AcousticPhonons(self.E)
		self.OP_a = super(g_Valley, self).OpticalPhonons(self.E,"A")
		self.OP_e = super(g_Valley, self).OpticalPhonons(self.E,"E")
		self.CS_a = super(g_Valley, self).CentralToSatellite(self.E,"A")
		self.CS_e = super(g_Valley, self).CentralToSatellite(self.E,"E")


class l_Valley(Valley):

	def __init__(self, E):

		super(l_Valley,self).__init__("L") 
		self.E = E
		self.AC = super(l_Valley,self).AcousticPhonons(self.E)
		self.OP_a = super(l_Valley,self).OpticalPhonons(self.E,"A")
		self.OP_e = super(l_Valley,self).OpticalPhonons(self.E,"E")
		self.SC_a = super(l_Valley,self).SatelliteToCentral(self.E,"A")
		self.SC_e = super(l_Valley,self).SatelliteToCentral(self.E,"E")
		self.SS_a = super(l_Valley,self).SatelliteToSatellite(self.E,"A")
		self.SS_e = super(l_Valley,self).SatelliteToSatellite(self.E,"E")


###################################################
# Elecrton Superclass
#
###################################################
class kVector(object):

	def __init__(self, kz, kr):
	
		self.kz  = kz 
		self.kr  = kr 

		self.mag = np.sqrt(kz**2+kr**2)

	def add(self, kVector):

		self.kz+=kVector.kz
		self.kr+=kVector.kr

		self.mag=np.sqrt(self.kz**2 + self.kr**2)


class Electron(GaAs, Universe):

	# Initialization method
	def __init__(self):

		GaAs.__init__(self) ; Universe.__init__(self);
		
		self._RANDOM_ = r.SystemRandom()

	# Methods to return the total momentum, the energy and the effective mass.
	def _getK(self,_E, _V): 
		return np.sqrt((2*self._getM(_V)*_E)/(self.hbar**2))

	def _getE(self,_K, _V): 
		return ((self.hbar**2)*(_K**2))/(2*self._getM(_V))

	def _getM(self,_V):

		if _V == "G": 
			return self.mc

		if _V == "L": 
			return self.ms

	# A method to increment the k-vector for free flight
	def _kField(self, ki, kf):

		return kVector( ki.kz + kf.kz , ki.kr + kf.kr)

	# Isotropic scattering events (component modification). To this
	# function we pass the FINAL ENERGY after scattering !
	def _kIsotropic(self, Ef, V):
		r1 = self._RANDOM_.random()
		kz = self._getK(Ef,V)*(1- (2*r1))
		kr = self._getK(Ef,V)*np.sqrt(4*r1*(1-r1))
		return kVector(kz, kr)

	# Anisotropic scattering events (component modification)
	def _kAnisotropic(self, Ef, Ei, V, ki):

		# Generate random numbers and calculate Xi
		r1 = self._RANDOM_.random()
		Xi = 2*np.sqrt(Ei*Ef)/(np.sqrt(Ei)-np.sqrt(Ef))**2

		# Calculate the various sines and cosines
		# for coordinate transformation
		cos_t = ((1+Xi)-np.power((1+2*Xi),r1))/Xi
		sin_t = np.sqrt(1-cos_t**2)
		cos_a = ki.kz/ki.mag
		sin_a = np.sqrt(1-cos_a**2)
		cos_p = math.cos(2*self.pi*r1)

		# Calculate the new k componenets
		kz = self._getK(Ef,V)*(-1*sin_a*sin_t*cos_p + cos_a*cos_t)
		kr = self._getK(Ef,V)*np.sqrt(1 - (-1*sin_a*sin_t*cos_p + cos_a*cos_t)**2)
		return kVector(kz, kr)


class _Electron(Electron):

	def __init__(self, _E, _V):
		
		super(_Electron, self).__init__()
		
		# Now we need to set the initial conditions. For the k-vector we calculate
		# a magnitude corresponding to the initial energy. Then we generate random
		# componenets by calling isotropic scattering
		self.E,self.T,self.V = [_E], [0.0], [_V]
		self.k = [super(_Electron,self)._kIsotropic(_E,_V)]

	# Methods to get the current state of the system
	def getE(self): 
		return self.E[-1]
	
	def getV(self): 
		return self.V[-1]
	
	def getT(self): 
		return self.T[-1]
	
	def getK(self): 
		return self.k[-1]

	# Methods to call the superclass methods on the current energy
	def calcE(self, _K): 
		return super(_Electron, self)._getE(_K, self.V[-1])
	
	def calcK(self, _E): 
		return super(_Electron, self)._getE(_E, self.V[-1])

	# Setter methods: These will append the corresponding value to vector.
	def addV(self,_V): 
		self.V.append(_V)
	
	def addE(self,_E): 
		self.E.append(self.E[-1] + _E)
	
	def addT(self,_T): 
		self.T.append(self.T[-1] + _T)

	def addK(self, mode, DELTA=None):

		if mode == "_Self": 
			self.k.append(self.k[-1])
		
		if mode == "Field": 
			self.k.append(super(_Electron,self)._kField(DELTA, self.k[-1]))
		
		if mode == "__Iso": 
			self.k.append(super(_Electron,self)._kIsotropic(self.E[-1], self.V[-1]))
		
		if mode == "anIso": 
			self.k.append(super(_Electron,self)._kAnisotropic(self.E[-1], self.E[-2], self.V[-1],self.k[-2]))


	# A method to return things out of the k-vector
	def returnK(self, mode):
		
		if mode is "mag": 
			return [i.mag for i in self.k]
		
		if mode is "kz": 
			return [i.kz for i in self.k]
		
		if mode is "kr": 
			return [i.kr for i in self.k]

	# A method to put everything together once the simulation finishes
	def buildData(self):
		
		self.xV = self.V[2::2]
		self.xM = [self._getM(i) for i in self.V[2::2]]
		self.xT = [self.T[i+1]-self.T[i] for i in range(len(self.T)-1)][0::2]
		self.VELOCITY = [self.hbar*kz/m for (kz,m) in zip(self.returnK("kz"),self.xM)]
		self.VMEAN = sum(self.VELOCITY)/len(self.VELOCITY)


# This is a class which only stores and calculates averages
class aElectron(Electron):

	def __init__(self, _E, _V):

		super(aElectron, self).__init__()

		self.E, self.T, self.V, self.counter = _E, 0.0, _V, 1
		self.k = super(aElectron,self)._kIsotropic(_E,_V)

		# In this class the velocity must be calculated explicitly
		# after each real scattering event
		self.v = self.hbar*self.k.kz/super(aElectron,self)._getM(self.V)

		# Time values to keep track of valley occupancy
		self.g, self.l = 0.0,0.0

	# Methods to call the superclass methods on the current energy
	def calcE(self, _K): 
		return super(aElectron, self)._getE(_K, self.V)
	
	def calcK(self, _E): 
		return super(aElectron, self)._getE(_E, self.V)

	# Get methods to return the current state of the electron
	def getV(self): 
		return self.V
	
	def getE(self): 
		return self.E
	
	def getT(self): 
		return self.T
	
	def getK(self): 
		return self.k

	# Set methods to increment the current state of the electron
	def addV(self, _V): 
		self.V =_V

	def addE(self, _E): 
		self.E+=_E
	
	def addT(self, _T):
		if self.V == "G": self.g+=_T; self.T+=_T
		if self.V == "L": self.l+=_T; self.T+=_T

	def addK(self, mode, DELTA = None):

		if mode == "_Self": 
			pass
		
		if mode == "Field":
			self.k = (super(aElectron,self)._kField(DELTA, self.k))
			self.v += self.hbar*self.k.kz/super(aElectron,self)._getM(self.V)
			self.counter+=1

		if mode == "__Iso":
			self.k = (super(aElectron,self)._kIsotropic(self.E, self.V))
			self.v += self.hbar*self.k.kz/super(aElectron,self)._getM(self.V)
			self.counter+=1
	
		if mode == "anIso":
			self.k = (super(aElectron,self)._kAnisotropic(self.E, self.E-DELTA, self.V,self.k))
			self.v += self.hbar*self.k.kz/super(aElectron,self)._getM(self.V)
			self.counter+=1

	def buildData(self):
		self.VMEAN = self.v/float(self.counter)
		self.VALLEY = {"G": self.g/self.T, "L": self.l/self.T}




###################################################
# The MONTE CARLO simulator
#
###################################################
class MonteCarlo(object):

	def __init__(self, _E, _FIELD, _nsteps, mode="FULL", U=Universe()):

		# Initial Conditions and Valleys
		self.E = _E
		self.n = int(_nsteps)
		self.FIELD = _FIELD
		self.counter = 0


		# Initialize the random number generator
		self._RANDOM_ = r.SystemRandom()

		# Initialize the mode (AVG or FULL)
		self.mode = mode

		# Sum over the Gammas to get the total Gamma. Also we
		# need the maximum gamma for the calculation of "self
		# scattering". These are NULL scattering events which
		# do not change the energy of the electron.
		G = g_Valley(self.E)
		L = l_Valley(self.E)

		GAMMA_G = G.AC + G.OP_e + G.OP_a + G.CS_e + G.CS_a
		GAMMA_L = L.AC + L.OP_e + L.OP_a + L.SC_e + L.SC_a + L.SS_e + L.SS_a

		# Total Gamma for each valley. This is needed for calculation
		# of self scattering and free flight times.
		self.GAMMA_G_max = max(GAMMA_G)
		self.GAMMA_L_max = max(GAMMA_L)

		# Normalization of the total Gamma Vector
		GAMMA_G = GAMMA_G/max(GAMMA_G)
		GAMMA_L = GAMMA_L/max(GAMMA_L)

		# Calculate the self Scattering Vectors
		self.GAMMA_Gss = [max(GAMMA_G)]*len(GAMMA_G) - GAMMA_G
		self.GAMMA_Lss = [max(GAMMA_L)]*len(GAMMA_L) - GAMMA_L

		# Store the Normalized Scattering Rates (Gamma-Valley)
		self.GAC = G.AC/self.GAMMA_G_max
		self.GOP_e = G.OP_e/self.GAMMA_G_max
		self.GOP_a = G.OP_a/self.GAMMA_G_max
		self.GCS_e = G.CS_e/self.GAMMA_G_max
		self.GCS_a = G.CS_a/self.GAMMA_G_max
		self.GSS = self.GAMMA_Gss

		# Store the Normalized Scattering Rates (L-Valley)
		self.LAC = L.AC/self.GAMMA_L_max
		self.LOP_e = L.OP_e/self.GAMMA_L_max
		self.LOP_a = L.OP_a/self.GAMMA_L_max
		self.LSC_e = L.SC_e/self.GAMMA_L_max
		self.LSC_a = L.SC_a/self.GAMMA_L_max
		self.LSS_e = L.SS_e/self.GAMMA_L_max
		self.LSS_a = L.SS_a/self.GAMMA_L_max
		self.LSS = self.GAMMA_Lss

		# Define the scattering matricies. These contains the probabilities
		# of all the various scattering events at each energy for both valleys
		# Given the energy, and the valley one is living in, one will generate
		# A scattering event. Each COLUMN represents an energy range and is
		# accessed via [:,i]. The sum of each column is unity.
		self.G_SCATTERING = np.asarray( [
			self.GSS, 
			self.GAC,
			self.GOP_e, 
			self.GOP_a, 
			self.GCS_e, 
			self.GCS_a
		])

		self.L_SCATTERING = np.asarray( [
			self.LSS,
			self.LOP_e,
			self.LSC_e,
			self.LSS_e,
			self.LAC,
			self.LOP_a,
			self.LSC_a,
			self.LSS_a
		])

		# Component vectors. This tells which way to calculate new components
		# given some scattering mechanism.
		self.G_COMP = np.asarray([
			"_Self", 
			"__Iso",
			"anIso", 
			"anIso",
			"__Iso", 
			"__Iso"
		])
		
		self.L_COMP = np.asarray([
			"_Self",
			"anIso",
			"__Iso",
			"__Iso",
			"__Iso",
			"anIso",
			"__Iso",
			"__Iso"
		])


		# An increment/decrement vector for the energy of the final state
		# after scattering.
		self.G_ENERGY = np.asarray([
			0.0, 
			0.0,
			-G.wOP*U.hbar, 
			G.wOP*U.hbar,
			-G.wE*U.hbar-G.D, 
			G.wE*U.hbar-G.D
		])

		self.L_ENERGY = np.asarray([
			0.0, 
			0.0,
			-G.wOP*U.hbar, 
			G.wOP*U.hbar,
			-G.wE*U.hbar+G.D, 
			G.wE*U.hbar+G.D,
			-G.wE*U.hbar, 
			G.wE*U.hbar
		])

		# Valley switching vectors. These contains the final state
		# of each scattering event (valley).
		self.G_VALLEY = np.asarray(["G","G","G","G","L","L"])
		self.L_VALLEY = np.asarray(["L","L","L","L","G","G","L","L"])

		# Initialize the electron object. This will set the initial energy E,
		# valley "G", k-vector and its components, and time 0. The calculation
		# of the initial components of k is obtained via the isotropic scattering.
		# formula.
		if self.mode == "FULL": 
			self.e = _Electron(0.01, "G")
		
		if self.mode == "AVG" : 
			self.e = aElectron(0.01, "G")

	def run(self, U=Universe()):

		while True:

			# Calculate the free flight time and kvector
			_T = self.GetTime()
			_K = kVector(-1*self.FIELD*_T/U.hbar, 0.0)

			# Update the state of the electron after free flight
			self.e.addT(_T)
			self.e.addV(self.e.getV())
			self.e.addK("Field", _K)
			self.e.addE(self.e.calcE(self.e.getK().mag)-self.e.getE())
	
			# Perform scattering Event
			self.Scatter()
			
			# If the number of scattering events has been reached
			# then break the simulation and return the electron
			if self.counter > self.n: 
				break
			
				return self.e
		
	# Determination of the free flight time
	def GetTime(self):
		
		if self.e.V[-1] == "G":
			return (-1.0/self.GAMMA_G_max)*np.log(self._RANDOM_.random())

		if self.e.V[-1] == "L":	
			return (-1.0/self.GAMMA_L_max)*np.log(self._RANDOM_.random())


	# Determination the the scattering Event
	def Scatter(self):

		# Get the index of the nearest energy to the energy vector we are considering
		_Ei = min(enumerate(self.E), key=lambda e: abs(e[1]-self.e.getE() ))[0]

		# Get the column corresponding to this energy in the scattering matrix
		# Note that this depends on the valley we happen to be in.
		if self.e.getV() == "G": 
			_S = self.G_SCATTERING[:,_Ei]
		
		if self.e.getV() == "L": 
			_S = self.L_SCATTERING[:,_Ei]
	
		_S = [sum(_S[0:i]) for i in range(len(_S)+1)]


		# Throw a random number to determine the scattering event.
		while True:
	
			num = r.random()
	
			for i,v in enumerate(_S):
	
				if _S[i]<num<_S[i+1]: 
					break

				if self.e.getV() == "G" and self.e.getE()+self.G_ENERGY[i]>0:
					break
				
				if self.e.getV() == "L" and self.e.getE()+self.L_ENERGY[i]>0:
					break

		# If the scattering is real increment counter
		if i!=0: 
			self.counter+=1 
			#print self.counter

		# Update the electron energy and valley
		if self.e.getV() == "G":
			#print self.e.getE()+self.G_ENERGY[i],self.G_ENERGY[i], self.e.getV(),i,self.G_COMP[i]
			self.e.addV(self.G_VALLEY[i])
			self.e.addT(0.0)
			self.e.addE(self.G_ENERGY[i])
			self.e.addK(self.G_COMP[i], self.G_ENERGY[i])
			return None
		
		if(self.e.getV()=="L"):
			#print self.e.getE()+self.L_ENERGY[i],self.L_ENERGY[i], self.e.getV(),i,self.L_COMP[i]
			self.e.addV(self.L_VALLEY[i])
			self.e.addT(0.0)
			self.e.addE(self.L_ENERGY[i])
			self.e.addK(self.L_COMP[i], self.L_ENERGY[i])
			return None


if __name__ == "__main__":
	
	##########################################
	#TESTING OF ELECTRON CLASS
	#
	##########################################
	if False:
	
		E, F ,NPOINTS = np.linspace(0.0,1.0,100), 1000, 10
		M = MonteCarlo(E, F, NPOINTS, "FULL")
		Electron = M.run()
		Electron.buildData()

		#plt.plot(Electron.T, Electron.returnK("mag"))
		plt.plot(Electron.T, Electron.E)
		plt.xlabel("Time $(s)$")
		plt.ylabel("Energy $(eV)$")
		plt.show()

	#####################################
	# PLOTS OF THE SCATTERING LIFETIMES #
	#####################################
	if True:
		E = np.linspace(0.00,1.0,100)
		G,L = g_Valley(E), l_Valley(E)
		
		# Acoustic Phonons
		fig = plt.figure(1)
		ax1 = fig.add_subplot(111)
		h1=ax1.semilogy(E,G.AC,"k")
		h2=ax1.semilogy(E,L.AC,"r")
		ax1.set_ylim(5e10,2e13)
		ax1.set_xlabel("Energy $(eV)$")
		ax1.set_ylabel("Acoustic Phonons $(s^{-1})$")
		plt.legend([h1,h2],["$\Gamma$ Valley","L valley"])
		plt.grid(which="both")
		
		# G Valley (2 plots)
		fig = plt.figure(2)
		ax1 = fig.add_subplot(111)
		h1=ax1.plot(E,G.OP_a,"k")
		h2=ax1.plot(E,G.OP_e,"r")
		ax1.set_ylim(0,1e13)
		ax1.set_xlabel("Energy $(eV)$")
		ax1.set_ylabel("$\Gamma$ Optical Phonons $(s^{-1})$")
		plt.legend([h1,h2],["Absorption", "Emission"])
		plt.grid(which="both")
		
		fig = plt.figure(3)
		ax1 = fig.add_subplot(111)
		h1=ax1.semilogy(E,G.CS_a,"k")
		h2=ax1.semilogy(E,G.CS_e,"r")
		ax1.set_ylim(1e12,3e14)
		ax1.set_xlabel("Energy $(eV)$")
		ax1.set_ylabel("$\Gamma$ Central-Satellite $(s^{-1})$")
		plt.legend([h1,h2],["Absorption", "Emission"])
		plt.grid(which="both")
		
		# L Valley (3 plots)
		fig = plt.figure(4)
		ax1 = fig.add_subplot(111)
		h1=ax1.plot(E,L.OP_a,"k")
		h2=ax1.plot(E,L.OP_e,"r")
		ax1.set_ylim(0,1.6e13)
		ax1.set_xlabel("Energy $(eV)$")
		ax1.set_ylabel("L Optical Phonons $(s^{-1})$")
		plt.legend([h1,h2],["Absorption", "Emission"])
		plt.grid(which="both")
		
		fig = plt.figure(5)
		ax1 = fig.add_subplot(111)
		h1=ax1.semilogy(E,L.SC_a,"k")
		h2=ax1.semilogy(E,L.SC_e,"r")
		ax1.set_ylim(1e11,1e13)
		ax1.set_xlabel("Energy $(eV)$")
		ax1.set_ylabel("L Satellite-Central $(s^{-1})$")
		plt.legend([h1,h2],["Absorption", "Emission"])
		plt.grid(which="both")
		
		fig = plt.figure(6)
		ax1 = fig.add_subplot(111)
		h1=ax1.semilogy(E,L.SS_a,"k")
		h2=ax1.semilogy(E,L.SS_e,"r")
		ax1.set_ylim(1e12,1e14)
		ax1.set_xlabel("Energy $(eV)$")
		ax1.set_ylabel("L Satellite-Satellite $(s^{-1})$")
		plt.legend([h1,h2],["Absorption", "Emission"])
		plt.grid(which="both")
		
		plt.show()
