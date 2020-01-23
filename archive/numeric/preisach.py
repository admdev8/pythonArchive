#!/usr/bin/env python 
import numpy as np 
import matplotlib.pyplot as plt

# Implementation of preisach kernel (archived)
class Preisach:

	def __init__(self, npoints = 100):

		pass

	# Method to expand domain
	def domain(self, domain):	

		return np.array( list(domain) + list(domain[::-1]) )
	
	# Method to evaluate an individual relay
	def relay(self, domain, r, d):

		# Evaluate ascending and descending branches
		asc = [-1. if _d <= ( r - d ) else 1. for _d in domain]
		des = [-1. if _d <= ( r + d ) else 1. for _d in domain[::-1]]

		# Return concatenated array
		return np.array( asc + des )

	# Method to evaluate a preisach kernel
	def kernel(self, domain, rmin, rmax, dmin, dmax, npoints = 100):

		# Kernel object prepare zeros
		kernel = np.zeros_like( self.domain(domain) )

		# Numeric range of relay positions
		for r in np.linspace(rmin, rmax, npoints):

			# Numeric range of relay widths
			for d in np.linspace(dmin, dmax, npoints):

				kernel += self.relay(domain, r, d)

		# Return kernel object
		return kernel / npoints**2


# Test the kernel
if __name__ == "__main__":

	# Create domain
	domain = np.linspace(0,1,1000)
	
	# Initialize class
	P = Preisach()

	_domain = P.domain(domain)
	for npoints in [2, 4, 8, 16, 32, 64]:

		_kernel = P.kernel(domain, 0.25, 0.75, 0.0, 0.2, npoints)
		plt.plot(_domain, _kernel)
	
	plt.show()
