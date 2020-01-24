#!/usr/bin/env python 
import sys
import collections

# Numerics
import cmath
import numpy as np

# Import utility functions
from ..utilites import *

# Generate complex number (utility)
def d2r(mag, ang): 
	return cmath.rect(mag, np.pi*ang/180)

# s2p data (network analyzer) (archive)
class s2p:

	def __init__(self,_dat, size): 

		# Store Size
		self.data = _dat
		self.size = 1

	@classmethod
	def fromfile(cls,fname):

		with open(fname) as f:
		content = f.readlines()
	
		_dat = collections.OrderedDict([('freq', []),('s11', []),('s12',[]),('s21', []),('s22',[])])
		for line in content:

			## Ignore comment section 
			if line[0] == '#' or line[0] == '!': 
				continue 

			_dat['freq'].append(float(line.split()[0]))
			_dat['s11'].append(complex(float(line.split()[1]),float(line.split()[2])))
			_dat['s21'].append(complex(float(line.split()[3]),float(line.split()[4])))
			_dat['s12'].append(complex(float(line.split()[5]),float(line.split()[6])))
			_dat['s22'].append(complex(float(line.split()[7]),float(line.split()[8])))

		size = len(_dat['freq'])
		d = cls(_dat, size) 
		return d


# s4p data (network analyzer) (archive)
class s4p:

	def __init__(self,_dat, _mod, size): 

		# Store Model
		self.data	= _dat
		self.model = _mod 

		# Store Size
		self.size	= size

		# Call other methods on data 
		self.data['K'], self.data['MSG'], self.data['MAG'], self.data['U'] = self.active_params(_dat)
		self.model['K'], self.model['MSG'], self.model['MAG'], self.model['U'] = self.active_params(_mod)


	@classmethod
	def fromfile(cls,fname):

		with open(fname) as f:
		content = f.readlines()
	
		_dat = collections.OrderedDict([('freq', []),('s11', []),('s12',[]),('s21', []),('s22',[])])
		_mod = collections.OrderedDict([('freq', []),('s11', []),('s12',[]),('s21', []),('s22',[])])
		 

		for i,line in enumerate(content):

			## Ignore comment section 
			if line[0] == '#' or line[0] == '!': 
				continue 

			# Ignore empty lines
			if len(line.split()) == 0: 
				continue

			# Data lines
			if len(line.split()) == 9: 

				_dat['freq'].append(float(line.split()[0]))
				_mod['freq'].append(float(line.split()[0]))

				## Load the data
				_dat['s11'].append(d2r(float(content[i].split()[1:][0]),float(content[i].split()[1:][1])))
				_dat['s12'].append(d2r(float(content[i].split()[1:][2]),float(content[i].split()[1:][3])))
				_dat['s21'].append(d2r(float(content[i+1].split()[0]),float(content[i+1].split()[1])))
				_dat['s22'].append(d2r(float(content[i+1].split()[2]),float(content[i+1].split()[3])))
								 

				 ## Load the model
				_mod['s11'].append(d2r(float(content[i+2].split()[4]),float(content[i+2].split()[5])))
				_mod['s12'].append(d2r(float(content[i+2].split()[6]),float(content[i+2].split()[7])))
				_mod['s21'].append(d2r(float(content[i+3].split()[4]),float(content[i+3].split()[5])))
				_mod['s22'].append(d2r(float(content[i+3].split()[6]),float(content[i+3].split()[7])))
								 

		size = len(_dat['freq'])
		d = cls(_dat,_mod, size) 
		return d

	# Method to calculate active circuit parameters (archive)
	def active_params(self,_dat):
		
		K, MSG, MAG, U = [],[],[],[]

		for i in range(self.size): 
		
			delta = (_dat['s11'][i]*_dat['s22'][i] - _dat['s12'][i]*_dat['s21'][i])
			K.append((1 - (abs(_dat['s11'][i]))**2 - (abs(_dat['s22'][i]))**2 + abs(delta)**2)/(2*abs(_dat['s12'][i]*_dat['s21'][i])));
			U.append(abs(_dat['s21'][i])**2/((1-(abs(_dat['s11'][i]))**2)*(1-(abs(_dat['s22'][i]))**2)))
			
			if K[-1] < 1:
				MSG.append(abs(_dat['s21'][i])/abs(_dat['s12'][i]));
				
			else: 
				MAG.append(abs(_dat['s21'][i])/abs(_dat['s12'][i])*(K[-1]-np.sqrt(K[-1]**2-1))) 

		MSG,MAG,U = [20*np.log10(i) for i in MSG],[20*np.log10(i) for i in MAG],[20*np.log10(i) for i in U]
		return K, MSG, MAG, U
