#!/usr/bin/python
import numpy as np

# Calculate the fermi level in graphene
def fermiLevel(T,carriersMeas=None): 

	k = 8.617e-5	   #eV/K
	hbar = 6.582e-16   #eV*s
	vfermi = 1e8	   #cm/s

	# Variable Constants
	eThermal = lambda T: k*T
	momentum = hbar*vfermi

	# Find carriers for various fermi levels
	carriers = []
	Ef = np.linspace(0.1,0.7,100)
	for i,E in enumerate(Ef):
		function = lambda r: r/(1+ np.exp(r-(E/eThermal(T))));
		n = (2/math.pi)*((eThermal(T)/momentum)**2)*integrate.quad(function,0,np.inf)[0]
		carriers.append(n)
  
	# Loop to find fermi level given carrier density
	if carriersMeas is not None:
		for i,n in enumerate(carriers):
			if carriersMeas < n: 
			   return Ef[i]

# Calculate the carriers given a fermiLevel
def fermiCarriers(T, Ef):
	
	# Define Constants
	k = 8.617e-5	   #eV/K
	hbar = 6.582e-16   #eV*s
	vfermi = 1e8	   #cm/s

	# Variable Constants
	eThermal = lambda T:k*T
	momentum = hbar*vfermi

	# Integrate over wavevectors
	function = lambda r: r/(1+np.exp(r-(Ef/eThermal(T))))
	carriers = lambda T: (2/math.pi)*((eThermal(T)/momentum)**2)*integrate.quad(function,0,np.inf)[0]

	# Return Carrier Density
	return carriers(T)
