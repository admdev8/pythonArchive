#!/usr/bin/env python
import math
import numpy as np
import numpy.linalg as la
from copy import deepcopy as c
import matplotlib.pyplot as plt

## For Plotting
from matplotlib import rc
rc("font",**{"family":"sans-serif","sans-serif":["Helvetica"],"size":20})
rc("text", usetex=True)
rc("xtick", labelsize=20)
rc("ytick", labelsize=20)
rc("legend", fontsize=20)

# WARNING: This script contains known errors (archive only)


# A container for the Silicon parameters.
class Si(object):

	def __init__(self):

		self.ni, self.ep, self.q, self.Vt = 1.5e10, 11.9*8.85e-12,1.602e-19,0.026

# A class to set up and solve Poisson via the finite difference method
class diodeSim(object):

	# To start we store the domain and range
	def __init__(self, d, scale=True):


		# Set up the domain and Range, as well as the differential
		# Operator Structure. The operator will later be modified
		# By the setting of boundary conditions via _embedBC.
		self.scale = np.sqrt( (0.01*Si.ep*Si.Vt)/(Si.q*Si.ni) )
		self.domain, self.size, self.h = np.asarray(d[1:-1])/self.scale,len(d)-2,(d[1]-d[0])/self.scale
		self.DD = self.ddOperator(self.size)


		# Set Flags for the boundary conditions. The solver will
		# Not run unless boundary conditions are set a BOTH ends.
		# We also need to add a BC vector for later.
		self.zth,self.nth, self.bcvector = False, False,[0]*(self.size)


		# Generate the vectors for the doping concentrations. These
		# will later be used for plotting
		self.Nd = np.asarray([self._Nd(i) for i in self.domain])
		self.Na = np.asarray([self._Na(i) for i in self.domain])
		
		# Define getter methods for the doping profile as a function
		# of position in the diode.
		def _Nd(self, d): return 1e18*np.exp(-(d*d*self.scale**2)/(1e-8))
		def _Na(self, d): return 1e15

		
		# Next we must consrtruct the discrete differential operator.
		# This is a tridiagonal matrix. This matrix is n+1 in size.
		# The zeroth and nth col/row must be modified to include the
		# boundary conditions !!
		def ddOperator(self, size):
	
			it, _dmat = 0,-2*np.eye(size)
	
			while it < size:
			
				try: 
			
					_dmat[it][it-1]=1
					_dmat[it][it+1]=1
		
				except IndexError:
					
					_dmat[0][-1]=0
					it+=1
			
			return (1.0/self.h**2)*(_dmat)
		
		# A method to embed the boundary conditions into the Operator
		# Matrix. This method MUST be called for the Newton simulator
		# to run.
		def _embedBC(self, _type, pos, value):
		if _type == "D":
		if pos == "0":
		self.bcvector[0] = (-1.0*float(value)/(self.h**2))
		self.zth = True
		if pos == "n":
		self.bcvector[-1] = (-1.0*float(value)/(self.h**2))


		self.nth = True
		if _type == "N":
		
			if pos == "0":
	
				self.DD[0][1] = (2.0/self.h**2)
				self.bcvector[0] = (-2.*float(value))/self.h
				self.zth = True
	
			if pos == "n":

				self.DD[-1][-2] = (2.0/self.h**2)
				self.bcvector[-1] = (-2.*float(value))/self.h
				self.nth = True

		# Method to return the "forcing term" at all points for a given
		# iteration step. Here s is a "Silicon" object. This allows for
		# the silicon parameters to be accounted for independently
		def getN(self, phi, s=Si()): 

			return np.asarray([np.exp(p) for p in phi])
	
		def getP(self, phi, s=Si()): 

			return np.asarray([np.exp(-p) for p in phi])
		
		def _Forcing(self, _r, s=Si()): 

			return (self.getN(_r) - self.getP(_r) - self.Nd/Si.ni + self.Na/Si.ni )

		def dForcing(self, _r, s=Si()): 

			return (self.getN(_r) + self.getP(_r))
		
		####################################
		# SOLVERS
		#
		def _NEWTON(self,guess, conv, omega=1.0):
	
			if self.zth and self.nth:
	
				while True:
	
				Jacobian, dF = c(self.DD), self.dForcing(guess)
	
				for i in range(self.size): 

					Jacobian[i][i]-=dF[i]
					F = np.dot(self.DD, guess) - self.bcvector - self._Forcing(guess)
					guess-=omega* np.dot(la.inv(Jacobian),F)
	

					#plt.plot(omega* np.dot(la.inv(Jacobian),F) )
					#print( max(abs(F)) )
				
					if max(abs(F)) < conv: 

						return guess

##############################
# SOLVING POISSON
#
if __name__ == "__main__":

	# First we need to set the domain vector and initialize the
	# constants for Silicon.
	Si, npoints = Si(), 300
	scale = (Si.ep*Si.Vt)/(Si.q*Si.ni)
	domain = [float(i) for i in np.linspace(0.0,10e-4,npoints)]
	_domain = [10000*i for i in domain[1:-1]]

	# Initialize the diodeSim object and then set the initial guess
	# Vector to something reasonable
	s = diodeSim(domain)
	f = (np.sqrt((s.Na - s.Nd)**2 + 4.*np.asarray([Si.ni]*len(s.Na))**2) + s.Nd - s.Na)/(2.0*Si.ni)
	guess, guess0 = [np.log(i) for i in f], [Si.Vt*np.log(i) for i in f]

	# Finally, embed he boundary conditions and run the simulation
	# ----> Mixed Conditions
	s._embedBC("D","0",guess[0]); s._embedBC("N","n",0.0)
	guess = Si.Vt*s._NEWTON(guess,1e-6)

	# ----> Two Dirchlet Conditions
	#s._embedBC("D","0",guess[0]); s._embedBC("D","n",guess[-1])
	#guess = Si.Vt*s._NEWTON(guess,1e-6)
	
	############################
	# PLOTTING
	#
	if False:
	
		fig = plt.figure(2)
		ax1 = fig.add_subplot(111)
		ax1.plot(_domain, guess,"r")
		ax1.plot(_domain, guess0,"r", linestyle="--")
		ax2 = ax1.twinx()
		ax2.semilogy(_domain, s.Nd,"k", linestyle="--")
		ax2.semilogy(_domain, s.Na,"k", linestyle="--")
		ax2.set_ylim(1e14,1e18)
		ax1.set_xlabel("Depth $(\mu m)$")
		ax1.set_ylabel("Potential $(V)$")
		ax2.set_ylabel("Carrier Density $(cm^{-3})$")
		ax1.set_ylim(-0.4,0.5)
		
		fig = plt.figure(3)
		ax1 = plt.subplot(111)
		ax2 = ax1.twinx()
		ax1.semilogy(_domain, Si.ni*s.getN(guess/Si.Vt),"k")
		ax2.semilogy(_domain, Si.ni*s.getP(guess/Si.Vt),"r")
		ax1.set_xlabel("Depth $(\mu m)$")
		ax1.set_ylabel("Electron Density $(cm^{-3})$")
		ax2.set_ylabel("Hole Density $(cm^{-3})$")
		plt.show()
