#!/usr/bin/env python 
import os

# Data parser utility functions
def drange(start, stop, step):
	
	r = start
	
	while r < stop:
		yield r
		r += step

def chunks(l, n):
	
	for i in xrange(0, len(l), n):

	yield l[i:i+n]

def nearest(array,value):
	
	idx = (np.abs(np.asarray(array)-value)).argmin()
	
	return array[idx],idx

def stream_lines(file_name):
	
	_file = open(file_name)
	
	while True:
	
	line = _file.readline()
	
	if not line:

		_file.close()
		break
	
	yield line

def gen_paths(_dir):

	for dirpath,_,filenames in os.walk(_dir):

		for f in filenames: 

			yield os.path.abspath(os.path.join(dirpath, f))
