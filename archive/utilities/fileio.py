#!/usr/bin/env python 
import os 

###########################
# 	FILE UTILITIES

# Method to check specific paths against what one is intersted in
def checkString(modifier, path): 
	
	for ck in modifier: 
		if ck in path:
			pass
		else:
			return None
 
	return True

# Search a directory for path(s) containing string
def find_path(root, modifier):
	
	for path, dirs, files in os.walk(root):
	
		for filename in files:
	
			if checkString(modifier, os.path.join(path, filename)):
	
				return os.path.join(path, filename)


# Generator function to stream lines out of a file
def stream_lines(file_name):
	
	_file = open(file_name)
	while True:
	
		line = _file.readline()
		if not line:
		
			_file.close()
			break
	
	yield line

# Read data from tab delimated file with single keyrow
# No header
def read_data(file_name):

	gen  = stream_lines(file_name)
	keys = gen.next().split()
	data = dict((k, []) for k in keys)

	while True:

		try:
			line = gen.next().split()
			
			for i,k in enumerate(keys):
				data[k].append(float(line[i]))
		
		except:
			return data
