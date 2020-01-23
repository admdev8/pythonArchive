#!/usr/bin/env python 
#
# This script contains the functions used for fitting velocity field data extracted 
# from pulsed-IV measurements in graphene microbrides. Scattering via surface optical 
# phonons is proposed as a dominant scattering mechaims
#
# A temperature dependent measurement of the carrier velocity vs. electric field 
# characteristic for as-grown and H-intercalated epitaxial graphene on SiC
#
# Journal of Applied Physics 113, 193708 (2013); https://doi.org/10.1063/1.4807162
#
import sys

# Scipy libs (for fitting)
import scipy.stats as stat
import scipy.optimize as optimize
import scipy.integrate as integrate

# Import ndfit (written)
import ndfit as ndf
import numpy as np



# This method is used to fit pulsed iv-data curves to tanh model 
def ivTanhFit(voltage,current,fit_guess,fit_step,convergence,maxdepth): 

	# Store Data
	_data = [voltage, current]

	# Hyperbolic model
	#fitfunc = lambda p,data: p[0] 
	fitfunc = lambda p,d: p[0]*(1 + p[1]*d[0])*np.tanh(p[2]*d[0]) 
	errfunc = lambda p,d: fitfunc(p,d)-d[1]
  
	# Run recursive fit with lattice throttling
	fit = ndf.ndfit(_data, fitfunc, errfunc, fit_guess, fit_step, convergence, maxdepth)
	fit.lattice_select("fco")
	fit.fit_throttling(1e-3, 2, 1e-3, 2)
	fit.fit_run()

	# Generate the nd_curve (the result)
	curve = fit.nd_curve()

	# Calculate Resistance
	rfunc = lambda p,v: p[0]*p[1]*(1+p[1]*v)+p[0]*p[2]*(1 + p[1]*v)*(1-np.tanh(p[2]*v)**2)  
	resistance = [ rfunc(fit.p_found, v) for v in curve[0]]
	resistance = np.divide(1,resistance) 

	# Print recursion informtiton 
	print "recursion depth:"+str(fit.depth) 
	print "convergence:"+str(fit.res_vec[-1])
	print "noise amplitude:"+str(fit.res_vec[-1]/len(voltage)) 
	print "params:"+str(fit.p_found) 
	return curve[0], curve[1], resistance, fit.p_found



# Fit velocity field data
def velocityFieldFit(field, velocity, fit_guess, fit_step, convergence, maxdepth, temp=None): 

	# Store Data 
	field = [float(e) for e in field]
	velocity = [float(v) for v in velocity]
	_data = [field[2::2], velocity[2::2]] 

	# Fitting to Velocity Field Model
	fitfunc = lambda p, d: (p[0]*d[0])/np.power((1+np.power((p[0]*d[0]/p[1]),p[2])), 1/p[2])					
	errfunc = lambda p, d: np.log10(fitfunc(p,d))-np.log10(d[1])					   
	
	# Run recursive fit
	fit = ndf.ndfit(_data, fitfunc, errfunc, fit_guess, fit_step, convergence, maxdepth)
	fit.lattice_select("fco")
	fit.fit_throttling(1e-1, 3, 1e-3, 3)
	fit.fit_run()

	curve = fit.nd_curve_extended(4)


	# Print recursion informtiton 
	print "recursion depth:"+str(fit.depth) 
	print "convergence:"+str(fit.res_vec[-1])
	print "noise amplitude:"+str(fit.noise_vec[-1]) 
	print "params:"+str(fit.p_found) 
	
	result = {}
	# Store Measured Data
	result['fieldMeas'] = field
	result['velocityMeas'] = velocity
	
	# Store Results of Fit
	result['fieldFit'] = curve[0]
	result['velocityFit'] = curve[1]
	result['mu'] = fit.p_found[0]
	result['vsat'] = fit.p_found[1]
	result['alpha'] = fit.p_found[2]

	# Return Result
	return result


# Fit surfact optical phonon energy for const carriers
def surfacePhononFit(carriers, temp, vsat, eFermi): 

	# Define Constants
	k = 8.617e-5	   #eV/K
	hbar = 6.582e-16   #eV*s
	vfermi = 1e8	   #cm/s
	momentum = hbar*vfermi
   
	# Functions for phonon calculation
	Noccup = lambda ePhonon,T: 1/(np.exp(ePhonon/eThermal(T))-1)
	const1 = lambda T: 2/(hbar*pi*sqrt(pi*_carriers(T))) 
	const2 = lambda T: 1/(4*pi*_carriers(T)*(momentum**2))
	
	# Get the carriers for a given Fermi energy
	eThermal = lambda T: k*T
	_carriers = lambda T: fermiCarriers(T,eFermi)

	# Calculate the Fermi Level
	eF_fit = eFermi

	############################################
	# Fit a Phonon Energy to Velocity vs. Temp #
	############################################
	_data = [temp,vsat] 

	fitfunc = lambda eP,d: const1(d[0])*eP[0]*sqrt(1-(eP[0]**2)*const2(d[0]))*(1/(Noccup(eP[0],d[0])+1))
	errfunc = lambda eP,d: fitfunc(eP,d) - d[1]
	fit_guess, fit_step, convergence, maxdepth = [0.100],[0.001],1e12,1000 
  
	# Run recursive fit
	fit = ndf.ndfit(_data, fitfunc, errfunc, fit_guess, fit_step, convergence, maxdepth)
	fit.fit_run()
	  
	# Print recursion informtiton 
	print ""
	print "--------- Phonon Fit ----------" 
	print "recursion depth:"+str(fit.depth) 
	print "convergence:"+str(fit.res_vec[-1])
	print "noise amplitude:"+str(fit.res_vec[-1]/len(temp)) 
	print "params:"+str(fit.p_found) 
 
	# Generate the nd_curve (the result)
	curve = fit.nd_curve_extended(2)

	# Get Best Fit value
	tempFit = curve[0]
	vsatFit = curve[1]
	
	result = {}   
	# Store Measured Data
	result['carriersMeas'] = carriers
	result['tempMeas'] = temp
	result['vsatMeas'] = vsat
	
	# Store Results of Fit
	result['tempFit'] = tempFit
	result['vsatFit'] = vsatFit
	result['ePhonon'] = fit.p_vec[-1]
	result['eFermi']  = eF_fit
	
	# Return Result
	return result
  

# Calculate carrier density for temperature and fermi energy
if __name__ == "__main__":
	
	print fermiCarriers(300, .250)
	print fermiCarriers(200, .250)
	print fermiCarriers(100, .250)
	print fermiCarriers(300, .520)
	print fermiCarriers(200, .520)
	print fermiCarriers(100, .520)