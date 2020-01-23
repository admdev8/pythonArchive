#!/usr/bin/env python 
import sys
import copy
import math
import operator
import itertools

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
	