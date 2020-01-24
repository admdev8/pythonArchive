#!/usr/bin/env python 
import re
import collections

# Numerics and plotting
import numpy as np

# Import utility functions
from ..utilites import *

## A class which holds two lists
class meas: 

	def __init__(self):

		## Make some emptylists
		self.v, self.i = [], []

	def appendv(self, _v, switch=1):

		self.v.append(switch*_v)
	
	def appendi(self, _i, switch=1):

		self.i.append(switch*_i)
	
# Measured transistor data from HP4145B semiconductor parameter analyzer (archive)
class transistorDC:
	
	def __init__(self,_dat, keys): 

		# Store Size
		self.data = _dat
		self.keys = keys

	@classmethod
	def VdId(cls,fname, switch=1):

		_dat = []
		with open(fname) as f:
			content = f.readlines()

		# Get out the control voltage. The re is needed to transform , into . 
		# for typecast of string literal to float
		size = len(content[0].split())/3
		_V1  = [ float(re.sub('[,]', '.', v)) for v in content[1].split()[0::3] ]

		# Create an OrderedDictered dict data structure. Each control voltage will 
		# have a meas object attached to it which will contain the 
		# corresponding IV curve
		_dat = collections.OrderedDict()
		for _v1 in _V1: 

			_dat[_v1] = meas()
		
		keys = _dat.keys()

		# Loop through the data
		for i,line in enumerate(content): 
			
			# A switch to prevent the header from being included
			if i == 0:

				continue 
			
			# Exclude empty lines
			if len(line.split()) == 0:

				continue

			# Data lines
			for col in range(size):

				# Fill up the meas objects with the data
				_dat[keys[col]].appendv(float(re.sub('[,]', '.',line.split()[3*col+1])),switch)
				_dat[keys[col]].appendi(float(re.sub('[,]', '.',line.split()[3*col+2]))/.1,switch)
				
		d = cls(_dat, keys) 
		return d

	@classmethod
	def VgId(cls,fname):

		_dat = []
		with open(fname) as f:
			content = f.readlines()

		# Get out the control voltage. The re is needed to transform , into . 
		# for typecast of string literal to float
		size = (len(content[0].split())-1)/2 
		_V1  = [ float(re.sub('[,]', '.', v)) for v in content[1].split()[1::2] ]

		# Create an OrderedDictered dict data structure. Each control voltage will 
		# have a meas object attached to it which will contain the 
		# corresponding IV curve
		_dat = collections.OrderedDict()
		
		for _v1 in _V1:

			_dat[_v1] = meas()
		
		keys = _dat.keys()

		# Loop through the data
		for i,line in enumerate(content): 
			
			# A switch to prevent the header from being included
			if i == 0: 
				continue 
			
			# Exclude empty lines
			if len(line.split()) == 0: 

				continue

			# Data lines
			for col in range(size):

				# Fill up the meas objects with the data
				_dat[keys[col]].appendv(float(re.sub('[,]', '.',line.split()[0])))
				_dat[keys[col]].appendi(float(re.sub('[,]', '.',line.split()[2*col+2]))/.1)

		d = cls(_dat, keys) 
		return d
