#!/usr/bin/env python 

# Import utility functions
from ..utilites import *

# Pulsed IV data (waveforms) 
class point:

	def __init__(self):

		self.vg  = 0.
		self.vd  = 0.
		self.WAV = {}

# Home built setup puled IV data (archive)
# Home built data format for waveforms
class pulsedIV:

	def __init__(self):

		pass

	def data(file_name, h=None): 
		
		## Check to see if we have a PulsedIV file
		gen = stream_lines(file_name)
		
		line = gen.next().split()
		if (line[0] != "##"): 
			return None
		
		line = gen.next().split()
		if (line[0] != "##"): 
			return None

		DAT = []
		gen = stream_lines(file_name)

		while True: 
			try:

				line = gen.next().split()
				if (line[0] == "##"):

					p = point()
					p.vg  = float(line[1])
					p.vd  = float(line[2])
				
					keys  = gen.next().split()[1:]
					p.WAV = dict((k, []) for k in keys) 

					while True:

						line = gen.next().split()
						if line[0] != "##": 

							[p.WAV[k].append(float(_)) for k,_ in zip(keys,line)]  

						else: 

							DAT.append(p); 
							break

			except StopIteration: return DAT


# Data from auriga pulsed IV setup
class auriga:

	def __init__(self):

		pass
	
	# Auriga pulsed IV data
	def data(file_name, compressed=False):

		gen  = stream_lines(file_name) 
		line = "" 
		h 	 = "!"

		if h is not None:

			while True:

				tmp  = line
				line = gen.next()
	
				if line.split() == [] or line.split()[0][0]==h: 

					continue

				else: 

					break
				
		# Have to take keys and first line as one thing due to the way 
		# Headers are arranged in Auriga files
		point = line.split(',')
		keys = [k.rstrip().replace(" ", "") for k in tmp.split(',')]

		# Create dictionary object (dict comprehension)
		data = dict((k, [float(p)]) for k,p in zip(keys, point))
		
		while True:
		
			try:

				line =  gen.next().split(",")

				for i,k in enumerate(keys):
	
					data[k].append(float(line[i]))

			except StopIteration: break
		
		if compressed: 

			return {"path": file_name, "v": data["VinNQ"], "i": data["IinNQ"]}
		else:
			
			return data
