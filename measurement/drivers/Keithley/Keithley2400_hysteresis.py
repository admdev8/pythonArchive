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

	def setVoltHx(self, _v1,_v2,_v3,_vcomp): 

		# Turn on Current Mode
		self.inst.write('SOUR:FUNC:MODE VOLT')

		# Set to Current Sweep Mode
		self.inst.write('SOUR:VOLT:MODE SWE')
		time.sleep(0.4)

		# Enable 4 wire sense 
		# self.inst.write('SYST:RSEN ON')

		# Enable voltage sense mode only and set compliance
		self.inst.write('SENS:FUNC "CURR:DC"')
		self.start, self.stop, self.step, self.comp = _v1, _v2, _v3, _vcomp

		self.inst.write('TRIG:COUN %s'%(int(((self.stop-self.start)/self.step))+1))
		time.sleep(0.4)
		self.inst.write('SOUR:DEL 0.5')
		time.sleep(0.4)
		self.inst.write('SENS:CURR:PROT %s'%(self.comp))
		
		# Enable some delay for the measurement to finish before 
		# reading out the data
		self.time = int((_v2-_v1)/_v3)

	def sweepUP(self):

		self.inst.write('SOUR:SWE:DIR UP')
		time.sleep(0.4)
		self.inst.write('SOUR:VOLT:STAR %s'%(self.start))
		time.sleep(0.4)
		self.inst.write('SOUR:VOLT:STOP %s'%(self.stop))
		time.sleep(0.4)
		self.inst.write('SOUR:VOLT:STEP %s'%(self.step))
	
	def sweepDOWN(self):

		self.inst.write('SOUR:SWE:DIR DOWN')
		time.sleep(0.4)
		self.inst.write('SOUR:VOLT:STAR %s'%(self.start))
		time.sleep(0.4)
		self.inst.write('SOUR:VOLT:STOP %s'%(self.stop))
		time.sleep(0.4)
		self.inst.write('SOUR:VOLT:STEP %s'%(self.step))

	def measureHalf(self):

		self.inst.write("OUTP ON")
		print "Stabilizing Current"
		time.sleep(2.0)
		dat = OrD({'Voltage':[],'Current':[] })

		# Start the hysteresis loop 
		self.sweepUP()
		self.inst.write("READ?")
		time.sleep(self.time)
		d = [float(i) for i in re.split('[,]',self.inst.read(8192))]

		dat['Voltage']+= d[0::5]
		dat['Current']+= d[1::5]
	   
		# Start the hysteresis loop 
		self.sweepDOWN()
		self.inst.write("READ?")
		time.sleep(self.time)
		d = [float(i) for i in re.split('[,]',self.inst.read(8192))]

		dat['Voltage']+= d[0::5]
		dat['Current']+= d[1::5]

		# Save the resistance calculations for last 
		dat['Resistance'] =  ndf.quotient(dat['Voltage'], dat['Current'])
		dat['dResistance'] =  ndf.derivative(dat['Voltage'], dat['Current'])

		# Turn off output and dump the data
		self.inst.write("OUTP OFF")
		pickle.dump(dat, open( self.fname, "wb" ) )

	def measureFull(self):

		self.inst.write("OUTP ON")
		print "Stabilizing Current"
		time.sleep(2.0)
		dat = OrD({'Voltage':[],'Current':[] })

		# Start the hysteresis loop 
		self.sweepUP()
		self.inst.write("READ?")
		time.sleep(self.time)
		d = [float(i) for i in re.split('[,]',self.inst.read(8192))]

		dat['Voltage']+= d[0::5]
		dat['Current']+= d[1::5]
	   
		# Start the hysteresis loop 
		self.sweepDOWN()
		self.inst.write("READ?")
		time.sleep(self.time)
		d = [float(i) for i in re.split('[,]',self.inst.read(8192))]

		dat['Voltage']+= d[0::5]
		dat['Current']+= d[1::5]

		# Invert Voltages 
		self.start, self.stop = -1*self.stop,-1*self.start

		# Start the hysteresis loop 
		self.sweepDOWN()
		self.inst.write("READ?")
		time.sleep(self.time)
		d = [float(i) for i in re.split('[,]',self.inst.read(8192))]

		dat['Voltage']+= d[0::5]
		dat['Current']+= d[1::5]
	   
		# Start the hysteresis loop 
		self.sweepUP()
		self.inst.write("READ?")
		time.sleep(self.time)
		d = [float(i) for i in re.split('[,]',self.inst.read(8192))]

		dat['Voltage']+= d[0::5]
		dat['Current']+= d[1::5]

		# Save the resistance calculations for last 
		dat['Resistance'] =  ndf.quotient(dat['Voltage'], dat['Current'])
		dat['dResistance'] =  ndf.derivative(dat['Voltage'], dat['Current'])

		# Turn off output and dump the data
		self.inst.write("OUTP OFF")
		pickle.dump(dat, open( self.fname, "wb" ) )


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
	fname = './data/dark/2x20/m2/graphene_hist_2000mV.pkl'

	# Initialize the Keithley
	k = Keithley2400(1, fname) 
	k.setVoltHx(0, 2.0, 0.025, 0.100)
	k.measureFull()

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
