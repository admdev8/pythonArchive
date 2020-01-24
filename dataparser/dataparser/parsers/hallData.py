#!/usr/bin/env python 
import re 

# Import utility functions
from ..utilites import *

# Parses Hall effect data (archive)
# Home built data format (archive)
class hallData: 
	
	def __init__(self, _dat, keys): 
		
		self.data = _dat
		self.keys = keys
		
	@classmethod 
	def fromfile(cls, fname): 
		
		with open(fname) as f: 

			content = f.readlines()

		regex = re.compile('\\*|\\#')

		# Get headers
		for line in content:

			if regex.findall(line.split()[0]) != []:

				continue

			else:

				keys = line.split()
				size = len(line.split())
				break

		# Establish dict with keys		
		_dat = OrD()
		for i in keys: 

			_dat[i] = [] 
		
		# Fill up the dict
		for line in content: 

			if line.split() == []: 

				continue
			
			if regex.findall(line.split()[0]) != []: 

				continue 
			
			for i,key in enumerate(_dat):

				if len(line.split()) != size: 

					continue
			   
				try:
					
					_dat[key].append(float(line.split()[i]))
				
				except:
				
					_dat[key].append(str(line.split()[i]))

		# Remove the key line from the data file
		for i,key in enumerate(_dat): 

			_dat[key] = _dat[key][1:]
		
		# Call class __init__ 
		d = cls(_dat, keys)
		return d
