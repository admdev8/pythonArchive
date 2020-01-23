#!/usr/bin/env python 

import numpy as np
from scipy import linspace,polyfit, polyval, poly1d, divide, mean

###########################
# 	NUMERIC OPERATIONS

# Method to calculate the taylor derivitive
def diffTaylor(v,i): 

	coeffs = polyfit(v,i,3)
	vFit = linspace(min(v),max(v),100) 
	iFit = polyval(coeffs, vFit) 
	polynomial = poly1d(coeffs)
   
	# Second order	
	didv = lambda ix:  3*coeffs[0]*((ix)**2) + 2*coeffs[1]*(ix) + coeffs[2]

	# Fifth order
	#didv = lambda ix: 
	#		6*coeffs[0]*((ix)**5) + 5*coeffs[1]*((ix)**4) + 4*coeffs[2]*((ix)**3) + 
	# 		3*coeffs[3]*((ix)**2) + 2*coeffs[4]*(ix) 	  + coeffs[5]

	rFit = divide(1, [float(didv(val)) for val in vFit] )
	return vFit, iFit, rFit

# Calculate Correlation coefficient
def correlation(dat1,dat2): 

	top=mean(np.multiply((dat1-mean(dat1)),(dat2-mean(dat2))))
	bot=np.std(dat1)*np.std(dat2)
	return top/bot


# Golden ratio recursive definition
def phi(x, n = 0):
	if n == 30:
		return 1
	else:
		return x + 1. /phi(x, n+1) 


# Example of defining simple fit function and error function for 
# cureve fitting. IN this case we fit y = a/x.
#	p: parameters for optimizing 
# 	c: numerical constants 
class simple_fit:

	def __init__(self):
		pass

	def fitfunc(dat, p, c):
		return p[0] / dat[0]

	def errfunc(dat,p, c): 
		return fitfunc(dat, p, c) - dat[1]


###########################
# 	ARRAY OPERATIONS

## Useful generator function for splitting lists
def chunks(l, n):
	for i in xrange(0, len(l), n):

		yield l[i:i+n]

# Find nearest element to value in array
def find_nearest(array, value):
	
	idx = (np.abs(array-value)).argmin()
	
	return idx, array[idx]
