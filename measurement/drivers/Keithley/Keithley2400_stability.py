#!/usr/bin/env python 
import Gpib
import time

# Data serializtion 
import cPickle as pickle

# Plotting
import matplotlib.pyplot as plt

# Data structure
from collections import OrderedDict as OrD

#Connect to the Keithley
class Keithley2400: 

	def __init__(self, gpib, fname):
		
		# Save the file to be written 
		self.fname = fname

		# Initialize the device and restore defaults
		self.inst = Gpib.Gpib(0,gpib)
		self.inst.write("*RST")

		# Query IDN to make sure it works
		self.inst.write("*IDN?")
		print self.inst.read(100) 

		self.inst.write("SOUR:CLE:IMM")

	def setCurrMode(self, _i, _vcomp): 
		# Set to Current Mode
		self.inst.write("SOUR:FUNC:MODE CURR")
		
		# Set Desired Current
		self.inst.write("SOUR:CURR:IMM:AMPL %s"%(_i))
		
		# Enable all measurement functions 
		# self.inst.write('SENS:FUNC "VOLT","CURR","RES"')
  
		# Enable 4 wire sense 
		# self.inst.write('SYST:RSEN ON')

		# Enable voltage sense mode only
		self.inst.write('SENS:FUNC "VOLT"')
		self.inst.write('SENS:VOLT:PROT %s'%(_vcomp))

	def measure(self, _t, _num):

		self.inst.write("OUTP ON")
		print "Stabilizing Current"
		time.sleep(2.0)
		dat = OrD({'Time':[],'Voltage':[],'Current':[],'Resistance':[]})

		for n in range(int(_num)):

			 tic = datetime.datetime.now()						 
			 if n%10 == 1:
				print "------------------------------------"
				print "Dumping Data: %s(s) elapsed (%s)"%(dat["Time"][-1],n)
				print "------------------------------------"
				pickle.dump(dat, open( self.fname, "wb" ) )
				dumpToDat(self.fname)


			 d = [float(i) for i in re.split('[,]',self.inst.read(512))]

			 dat['Voltage'].append(d[0])
			 dat['Current'].append(d[1])
			 dat['Resistance'].append(dat['Voltage'][-1]/dat['Current'][-1])

			 time.sleep(_t)
			 print dat['Voltage'][-1], dat['Current'][-1], dat['Resistance'][-1]
			 
			 toc = datetime.datetime.now()
			 _time = toc-tic
			 _time.seconds + _time.microseconds/1e6
			 if dat["Time"] == []:
				 dat['Time'].append(_time.seconds+_time.microseconds/1e6)
			 else:
				 dat['Time'].append(dat["Time"][-1]+_time.seconds+_time.microseconds/1e6)
			 
 
		pickle.dump(dat, open( self.fname, "wb" ) )
		self.inst.write("OUTP OFF")

def dumpToDat(fname):

	# Load the .pkl file
	dat = pickle.load( open( fname, "rb" ) )
	dname = re.sub('.pkl','.dat', fname)
	f = open(dname, 'wb')

	# Write the keys to the first line
	keys = dat.keys()
	f.write('%s\t%s\t%s\t%s\n'%(str(keys[0]),str(keys[1]),
								str(keys[2]),str(keys[3])))

	# Write the data to later lines
	for i in range(len(dat[keys[0]])):
		f.write('%s\t%s\t%s\t%s\n'%(str(dat[keys[0]][i]),str(dat[keys[1]][i]),
									str(dat[keys[2]][i]),str(dat[keys[3]][i])))

	f.close()

if __name__ == "__main__": 

	# Set filename for data
	fname = './data/GaN/str1/80mA_24h.pkl'

	# Initialize the Keithley
	k = Keithley2400(1, fname) 
	k.setCurrMode(0.080, 20.0)
	k.measure(15, 6000)

	# Plot the data
	dat = pickle.load( open( fname, "rb" ) )
	plt.plot(dat['Time'], dat['Resistance'])
	plt.xlabel("Time (s)")
	plt.ylabel("Resistsnace (Ohm)")
	plt.savefig(re.sub('.pkl','.png', fname))
	plt.show()
	