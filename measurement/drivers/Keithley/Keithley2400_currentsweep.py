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

	def setCurrSweep(self, _i1,_i2,_i3,_vcomp): 

		# Turn on Current Mode
		self.inst.write('SOUR:FUNC:MODE CURR')

		# Set to Current Sweep Mode
		self.inst.write('SOUR:CURR:MODE SWE')
		time.sleep(0.4)
		self.inst.write('SOUR:CURR:STAR %s'%(_i1))
		time.sleep(0.4)
		self.inst.write('SOUR:CURR:STOP %s'%(_i2))
		time.sleep(0.4)
		self.inst.write('SOUR:CURR:STEP %s'%(_i3))
		time.sleep(0.4)
		self.inst.write('TRIG:COUN %s'%(int((_i2-_i1)/_i3)))
		time.sleep(0.4)
		self.inst.write('SOUR:DEL 0.5')

		# Enable 4 wire sense 
		# self.inst.write('SYST:RSEN ON')

		# Enable voltage sense mode only and set compliance
		self.inst.write('SENS:FUNC "VOLT"')
		self.inst.write('SENS:VOLT:PROT %s'%(_vcomp))
		
		self.time = int((_i2-_i1)/_i3)/1.4

	def measure(self):

		self.inst.write("OUTP ON")
		print "Stabilizing Current"
		time.sleep(2.0)
		dat = OrD()

		self.inst.write("READ?")
		time.sleep(self.time)

		d = [float(i) for i in re.split('[,]',self.inst.read(10000))]

		dat['Voltage'], dat['Current'] = d[0::5],d[1::5]
		dat['Resistance'] =  ndf.quotient(dat['Voltage'], dat['Current'])
		dat['dResistance'] =  ndf.derivative(dat['Voltage'], dat['Current'])

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
	fname = './data/GaN/str2/GaN_sweep_100mA.pkl'

	# Initialize the Keithley
	k = Keithley2400(1, fname) 
	k.setCurrSweep(0.001,0.100,0.001, 100.0)
	k.measure()

	# Dump to dat file
	dumpToDat(fname)

	# Plot the data
	dat = pickle.load( open( fname, "rb" ) )
	
	fig = plt.figure(1)
	ax1 = fig.add_subplot(111)
	ax1.plot(dat['Voltage'], dat['Current'])
	ax2 = ax1.twinx()
	ax2.plot(dat['Voltage'], ndf.derivative(dat['Voltage'], dat['Current']),'r-')
	ax2.plot(dat['Voltage'], ndf.quotient(dat['Voltage'], dat['Current']),'g-')

	ax1.set_xlabel("Voltage (V)")
	ax1.set_ylabel("Current (A)")
	ax2.set_ylabel("Resistance (Ohm)")
	plt.savefig(re.sub('.pkl','.png', fname))
	plt.show()

	




