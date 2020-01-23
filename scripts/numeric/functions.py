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

###########################
# 	SIGNAL OPERATIONS

# Method to despike a data set using convolution
def despike(counts, wavelist = [5,10,15,20] ): 

	# copy input arrry
	counts_despiked = copy.deepcopy(counts)

	# Loop through wavelist
	for val in waveList:

		# Create sinc function
		_wx = np.linspace(-1,1,val)
		_wy = 2*np.sinc(2*_wx)-np.sinc(_wx)

		# Perform convolution
		y = np.convolve(_wy/_wy.sum(), counts, mode='same')

		# Look at sign of convolution product with counts
		window = 3
		for i,val in enumerate( np.sign(counts*y) ): 

			# Negative values indicate a spike
			if val == -1:

				# Try to replace spike with surrounding data
				if i < len(counts_despiked) - window:
			
					counts_despiked[i] = (counts_despiked[i+window] + counts_despiked[i-window])/2


# Method to smooth data
def simple_smooth(self, data): 

	window_len, beta = 30, 30
	w = np.kaiser(window_len, beta)
	return np.convolve(w/w.sum(), data, mode='same')



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


# Linefit method (archived)
def linefit(x,y):

	_x = sum(x)/len(x)
	_y = sum(y)/len(y)

	dat,N,D = zip(x,y), 0.0, 0.0
	for pt in dat:
		N += (pt[0]-_x)*(pt[1]-_y) 
		D += (pt[0]-_x)*(pt[0]-_x) 

	return [N/D, _y - (N/D)*_x]


# Original nonlinear fit algorithm (archived)
def nonlinearFit(x,y,fitfunc,errfunc,fit_guess,fit_step, convergence=1e-6): 

	# Walk around parameter space
	res_vec, p_vec = [],[]

	#First need to set inital (p) and step size (d)
	p = [val for val in fit_guess]
	d = [val for val in fit_step]

	# Calculate the initial residual
	res = 0
	for n,val in enumerate(y): 
		res+= errfunc(p,x[n],y[n])**2
	res_vec.append(res)
	p_vec.append(p)

	c = 1
	nsteps = 10000
	while c < nsteps: 
  
		# 2n sides for an n-dimensional box in parameter space
		tmp_res,tmp_p = [],[]
		for i in range(2*len(p)): 

			# Positive Direction Steps
			if i<len(p):
				
				####### EDGES ############
				# Deep copy of p and index 
				pt = copy.deepcopy(p)
				pt[i]+=d[i]
				# Calculate residuals
				res = 0
				for n,val in enumerate(y):
					res+= errfunc(pt,x[n],y[n])**2
					
				# Save residual and p
				tmp_res.append(res)
				tmp_p.append(pt)
			
				####### CORNERS###########
				# Deep copy of p and index 
				pt = copy.deepcopy(p)
				pt = [n+d[i] for n in pt]
				# Calculate residuals
				res = 0
				for n,val in enumerate(y):
					res+= errfunc(pt,x[n],y[n])**2
				
				# Save residual and p
				tmp_res.append(res)
				tmp_p.append(pt)


			# Negative Direction Steps
			else:	  
			  
				####### EDGES ############
				# Deep copy of p and index 
				pt = copy.deepcopy(p)
				pt[i-len(p)]-=d[i-len(p)]
				# Calculate residuals
				res = 0
				for n,val in enumerate(y): 
					res+= errfunc(pt,x[n],y[n])**2
				
				# Save residual and p
				tmp_res.append(res)
				tmp_p.append(pt)
				
				####### CORNERS###########
				# Deep copy of p and index 
				pt = copy.deepcopy(p)
				pt = [n-d[i-len(p)] for n in pt]
				# Calculate residuals
				res = 0
				for n,val in enumerate(y): 
					res+= errfunc(pt,x[n],y[n])**2
				
				# Save residual and p
				tmp_res.append(res)
				tmp_p.append(pt)
								   
		# Save minimum values
		res_vec.append(min(tmp_res))
		p_vec.append(tmp_p[tmp_res.index(min(tmp_res))])
			 
		if (np.abs(res_vec[-1]) > np.abs(res_vec[-2])) and (min(res_vec)<convergence):
			p_found =  p_vec[-1]
			break

		else: 
			# Make the step
			p = p_vec[-1]
			c+=1

	if c == nsteps: 
		print "Warning: Probably not converged - try again with larger step"
		print "Returning final value" 
		p_found = p_vec[-1] 

	_x = np.linspace(min(x), max(x), 1000)
	_y = fitfunc(p_found,np.linspace(min(x), max(x), 1000))

	print p_found
	result = {}
	result["domain"] = _x
	result["range"] = _y 
	result["p_found"] = p_vec[-1]
	result["p_walk"] = p_vec
	result["residual"] = res_vec
	return result
	
	
###########################
# 	ARRAY OPERATIONS

## Useful generator function for splitting lists
def chunks(l, n):

	for i in xrange(0, len(l), n):

		yield l[i:i+n]

# Extract every other point
def zigzag(seq): 

	return seq[::2], seq[1::2]

# Find nearest element to value in array
def find_nearest(array, value):
	
	idx = (np.abs(array-value)).argmin()
	
	return idx, array[idx]


# Get exact value in array and index
def find_exact(self, array, v):

	try:
		
		ix = array.index(v)
		return (array[ix], ix)

	except ValueError:
	
		print("Value Not Found")
		return None
