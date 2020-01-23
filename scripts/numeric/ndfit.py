# nd_Fit.py version 0.0.2
#
# Welcome to the n-dimensional curve fit algorithm. This fits curves via 
# a least squares method. It can fit an m-dimensional curve to a fit 
# fit function of n-parameters. The m-dimensional curve must take the 
# form y = f(a1,...an : x1,...,xn). The algorithm requires the following:
# 
# 1) An initial guess for the parameters
# 2) Lattice constants of the fit lattice
# 3) A fit function and error function  
#
# The fit lattice is either face centered rhombohedral or body centered 
# rombohedral. Support for other fit lattice geometries is currently 
# under development.
#
# Author: M. Winters

#!/usr/bin/env python 
import sys
import copy
import math
import operator
import itertools

# Original ndfit fit routine as class (archived)
class ndfit(object): 
	
	def __init__(self, _dat, fitfunc, errfunc, _guess, _step, converge=1e-6, maxd=10000, l=100):

		#Set Recursion Parameters
		self.depth, self.maxdepth = 0,maxd
		sys.setrecursionlimit(maxd)

		# Save data convergence limit and points to generate for output
		self._dat, self.converge, self.l  = _dat, converge, l

		# Save error function and fit function 
		self.errfunc, self.fitfunc = errfunc, fitfunc

		#First need to set inital (p) and step size (d)
		self.p = [val for val in _guess]
		self.d = [val for val in _step]
			
		# And initialize lists for the recursive walk
		self.res_vec, self.noise_vec, self.p_vec = [],[],[]
		
		# Calculate Initial residual and noise
		self.p_vec.append(self.p)
		self.res_vec.append(self.nd_Residual(self.p))
		self.noise_vec.append(self.res_vec[-1]/len(_dat[0]))

	# A method to select lattices
	def lattice_select(self,name): 
		if name == "fco": 
			self.nd_fco_Lattice() 
		elif name == "bco": 
			self.nd_bco_Lattice() 
		else: 
			print "Illegal Lattice: Choose 'bco' or 'fco'"
			exit(1)

	# A method to set the parameters for Lattice throttle
	def fit_throttling(self, h_noise=1, h_factor=1, l_noise=1, l_factor=1):
		# Set Throttling Parameters
		self.accel = [h_noise,h_factor,l_noise,l_factor]

	# An internal method for lattice throttling
	def lattice_Throttle(self): 
		if self.noise_vec[-1] > self.accel[0]:
			return [map(operator.mul,i,[self.accel[1]]*len(self.p)) for i in self.lattice]
		elif self.noise_vec[-1] < self.accel[2]:
			return [map(operator.div,i,[self.accel[3]]*len(self.p)) for i in self.lattice]
		else: 
			return self.lattice
	
	# The run script itself
	def fit_run(self):

		# Set default (no lattice throttling) 
		if not hasattr(self, 'accel'):
			self.fit_throttling()

		# Build the default nd_Lattice (body-centered-orthorhombic)
		if not hasattr(self, 'lattice'): 
			self.nd_bco_Lattice() 

		# Pass to the recursive algorithm and do the fit
		self.nd_RecursiveFit(self._dat,self.fitfunc,self.errfunc,self.p,self.d,self.converge)

		# And save the found value for easy access
		self.p_found = self.p_vec[-1]

	# Generate standard output values of the fit
	def nd_curve(self): 

		# Initialize arrays for the domain 
		self.nd_curve = [map(lambda b:min(i)+(max(i)-min(i))*b/self.l,range(self.l)) for i in self._dat[:-1]]

		# Generate the range based on the domain as calculated above
		self.nd_curve.append([self.fitfunc(self.p_found,_tuple) for _tuple in zip(*self.nd_curve)])

		# return the nd_curve
		return self.nd_curve

	# Generate extended output values of the fit
	def nd_curve_extended(self, factor): 

		# Initialize arrays for the domain 
		self.nd_curve = [map(lambda b:(min(i)-min(i)/factor)+ factor*(max(i)-min(i))*b/self.l,range(self.l)) for i in self._dat[:-1]]

		# Generate the range based on the domain as calculated above
		self.nd_curve.append([self.fitfunc(self.p_found,_tuple) for _tuple in zip(*self.nd_curve)])

		# return the nd_curve
		return self.nd_curve

	# A method to generate an nd_residual with m variable and n parameters
	def nd_Residual(self,pt): 
		return sum([self.errfunc(pt,_tuple)**2 for _tuple in zip(*self._dat)])

	def nd_fco_Lattice(self): 
		
		self.lattice = []
		####### EDGES (2n) ##############
		nI = [[0]*len(self.p) for i in range(len(self.p))]
		for i in range(len(nI)): 
			nI[i][i]+=1
		_pfaces = [map(lambda a,b: a*b, vec,self.d) for vec in nI]

		nI = [[0]*len(self.p) for i in range(len(self.p))]
		for i in range(len(nI)): 
			nI[i][i]-=1
		_nfaces =  [map(lambda a,b: a*b, vec,self.d) for vec in nI]
		
		####### CORNERS (2^n) ###########
		# Need +++,-++,+-+,++-,+--,-+-,--+,---
		signList = list(itertools.product([1,-1],repeat=len(self.p)))
		_corners = [map(lambda a,b: a*b, vec,self.d) for vec in signList]
		self.lattice = _nfaces+_pfaces+_corners  
			
	def nd_bco_Lattice(self): 
		
		self.lattice = []
		####### CORNERS (2^n) ###########
		signList = list(itertools.product([1,-1],repeat=len(self.p)))
		self.lattice = [map(lambda a,b: a*b, vec,self.d) for vec in signList]

	# Take a walk around parameter space
	def nd_RecursiveFit(self,_dat,fitfunc,errfunc,fit_guess,fit_step,converge): 
					   
		tmp_res,tmp_p = [],[]
		# Calculate the residual for each lattice point
		for point in self.lattice_Throttle():
		  
			# Create point 
			pt = copy.deepcopy(fit_guess)
			pt = map(lambda a,b: a+b,pt,point)
		   
			# Save residual and p
			tmp_res.append(self.nd_Residual(pt))
			tmp_p.append(pt)
			
		# Save min residual and corresponding lattice point
		self.res_vec.append(min(tmp_res))
		self.noise_vec.append(math.sqrt(min(tmp_res))/len(_dat[0]))
		self.p_vec.append(tmp_p[tmp_res.index(min(tmp_res))])

		# Check the case of convergence
		if (self.res_vec[-1] > self.res_vec[-2]) and (min(self.res_vec)<converge):
			return 1

		# Check the case where maximum number of recursive calls is exceeded
		elif self.depth > self.maxdepth-100:
			print "Warning: Potentially not converged"
			return 2

		# Make the recursive call. This corresponds to taking a step in the 
		# parameter lattice towards a residual minimizing condition
		else: 
			self.depth +=1
			self.nd_RecursiveFit(_dat,fitfunc,errfunc,self.p_vec[-1],self.d, converge)	
