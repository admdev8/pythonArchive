#!/usr/bin/env python 
import Gpib
import time
import datetime
import cPickle as pickle

import re
import sys

import numpy as np
import matplotlib.pyplot as plt
 
#Connect to the Keithley
class Keithley_2400: 

	def __init__(self, gpib):
		
		# Initialize the device and restore defaults
		self.inst = Gpib.Gpib(0,gpib)
		self.inst.write("*RST")

		# Query IDN to make sure it works
		self.inst.write("*IDN?")
		print self.inst.read(100) 

		self.inst.write("SOUR:CLE:IMM")
		
	def setVoltMode(self, _v, _icomp): 
		# Set to Current Mode
		self.inst.write("SOUR:FUNC:MODE VOLT")
		
		# Set Desired Current
		self.inst.write("SOUR:VOLT:IMM:AMPL %s"%(_v))
		
		# Enable voltage sense mode only
		self.inst.write('SENS:FUNC "CURR"')
		self.inst.write('SENS:CURR:PROT %s'%(_icomp))
	   
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

	def outputsOn(self):
		self.inst.write("OUTP ON")
	def outputsOff(self):
		self.inst.write("OUTP OFF")


class Keithley_2000(object):
	def __init__(self,gpib):

		self.meter = Gpib.Gpib(0,gpib)
		self.meter.write("*RST")

		# Query IDN to make sure it works
		self.meter.write("*IDN?")
		print self.meter.read(100) 
		
	def _Imode(self): 
		self.meter.write("CONF:CURR:DC")
	
	def _Vmode(self): 
		self.meter.write("CONF:VOLT:DC")
	
	def getValue(self): 
		self.meter.write("READ?"); return float(self.meter.read(128))
 

# Driver Code Here
class source_66xx:

    def __init__(self,gpib):
        self.source = Gpib.Gpib(0,gpib)
      
        # Zero Voltages
        self.source.write("VSET 1,0")
        self.source.write("VSET 2,0")
        # Outputs off
        self.source.write("OUT 1,0")
        self.source.write("OUT 2,0")

    def setVoltage(self,v1set,v2set,compliance):
        if v1set > 40:
            print "V1 is too large"
        elif v2set > 40:
            print "V2 is too large"
        else:
            try:
                # Output One
                self.source.write("VSET 1,"+str(v1set))
                self.source.write("ISET 1,"+str(compliance))
                # Output Two
                self.source.write("VSET 2,"+str(v2set))
                self.source.write("ISET 2,"+str(compliance))
            except:
                print "Error: Check GPIB connection and address"

    def outputsOn(self):

        self.source.write("OUT 1,1")
        self.source.write("OUT 2,1")
    def outputsOff(self):
        
        self.source.write("OUT 1,0")
        self.source.write("OUT 2,0")
    
    # Getter methods
    def getVoltageCh1(self): 
        
        try: 
            self.source.write("VOUT? 1"); return float(self.source.read(64))
        
        except ValueError:
        
            print "Trying Again"
            self.source.write("VOUT? 1"); return float(self.source.read(64))

    def getVoltageCh2(self): 
        
        try: 
        
            self.source.write("VOUT? 2"); return float(self.source.read(64))
        
        except ValueError:
        
            print "Trying Again"
            self.source.write("VOUT? 2"); return float(self.source.read(64))


class Measurement:

	# Basic measurement class
	def __init__(self, _t, _num, _fname):
		
		# Save the measurement params
		self.t, self.num, self.fname = _t, _num, _fname
	  
	# Pass in the instrument handles. We need to do this so 
	# We can read out their data as we measure. 
	def measure(self,_Vgs, _Vds, _Ids):
		
		try: 
			if not isinstance(_Vgs,Keithley_2400): 
				raise ValueError

			if not isinstance(_Vds,source_66xx):
				raise ValueError

			if not isinstance(_Ids,Keithley_2000):
				raise ValueError

			self._Vgs, self._Vds, self._Ids = _Vgs,_Vds,_Ids
			self._meas()
		
		except ValueError:
			print "Insturments not set correctly: Try again"

	# The internal measurement code. This is a separate method to 
	# allow for type checking step. 
	def _meas(self):

		# Turn ON outputs
		_Vgs.outputsOn();_Vds.outputsOn()
		_VDS = self._Vds.getVoltageCh1()

		# Start the measurement sequence
		print "Stabilizing Bias"
		time.sleep(2.0)
		dat = OrD({'Time':[],'Vds':[],'Ids':[],'Rds':[],'Vgs':[],'Igs':[],'Rgs':[]})

		for n in range(int(self.num)):

			tic = datetime.datetime.now()						 
			if n%10 == 1:
				print "---------------------------------------------------------------"
				print "Dumping Data: %s(s) elapsed (%s)"%(dat["Time"][-1],n)
				print "---------------------------------------------------------------"
				pickle.dump(dat, open( self.fname, "wb" ) )
				#self.dumpToDat(self.fname)
				  
			# Get the Vds data
			d = [float(i) for i in re.split('[,]',self._Vgs.inst.read(512))]
			dat['Vgs'].append(d[0]); dat['Igs'].append(d[1]); dat['Rgs'].append(dat['Vgs'][-1]/dat['Igs'][-1])
			
			# Get the Vgs data
			dat['Vds'].append(_VDS); dat['Ids'].append(self._Ids.getValue())  
			dat['Rds'].append(dat['Vds'][-1]/dat['Ids'][-1])

			# Print the current data series
			print dat['Vds'][-1], dat['Ids'][-1], dat['Rds'][-1],dat['Vgs'][-1], dat['Igs'][-1], dat['Rgs'][-1]
			 
			# Sleep for Required Time
			time.sleep(self.t)
			toc = datetime.datetime.now()
			_time = toc-tic
			_time.seconds + _time.microseconds/1e6
			if dat["Time"] == []:
				dat['Time'].append(_time.seconds+_time.microseconds/1e6)
			else:
				dat['Time'].append(dat["Time"][-1]+_time.seconds+_time.microseconds/1e6)
			 
		pickle.dump(dat, open( self.fname, "wb" ))
		self.dumpToDat(self.fname)
		_Vgs.outputsOff()
		_Vds.outputsOff()

	def dumpToDat(self, fname):

		# Load the .pkl file
		dat = pickle.load( open( fname, "rb" ) )
		dname = re.sub('.pkl','.dat', fname)
		f = open(dname, 'wb')

	   # Write the keys to the first line
		keys = dat.keys()
		f.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(str(keys[0]),str(keys[1]),
												str(keys[2]),str(keys[3]),
												str(keys[4]),str(keys[5]),
												str(keys[6])))

		# Write the data to later lines
		for i in range(len(dat[keys[0]])):
			f.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(str(dat[keys[0]][i]),str(dat[keys[1]][i]),
													str(dat[keys[2]][i]),str(dat[keys[3]][i]),
													str(dat[keys[4]][i]),str(dat[keys[5]][i]),
													str(dat[keys[6]][i])))
		f.close()

if __name__ == "__main__": 

	# Set filename for data and initialize insturments
	_Vgs, _Vds, _Ids = Keithley_2400(1), source_66xx(2), Keithley_2000(3)

	Vglist,sgn,Vd = [0.,1.,2.,3.,4.,5.,6.,7.,8.], (-1,'m'), 1.
	
	for Vg in Vglist: 

		# Set Gate/Drain Voltage
		_Vgs.setVoltMode(sgn[0]*Vg,0.001)
		_Vds.setVoltage(Vd,0,0.100)
		_Ids._Imode()

		# Set up the measurement
		fname = './data/graphene/sampleA/device5/m2/%s%sV%sV.pkl'%(sgn[1],str(int(Vg)),str(int(Vd)))
		print 'Writing for %s%sV ;)'%(str(sgn[1]),str(Vg))
		print fname
		
		m = Measurement(0.0, 5000, fname)
		m.measure(_Vgs,_Vds, _Ids)
		
		#Set Outputs Off (Just in Case)
		_Vgs.outputsOff()
		_Vds.outputsOff()
		time.sleep(60)

	# Plot the data
	#dat = pickle.load( open( fname, "rb" ) )
	#plt.plot(dat['Time'], dat['Rds'])
	#plt.xlabel("Time (s)")
	#plt.ylabel("Resistsnace (Ohm)")
	#plt.savefig(re.sub('.pkl','.png', fname))
	#plt.show()