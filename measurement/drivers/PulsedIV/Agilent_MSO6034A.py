#!/usr/bin/env python 
#import pyvisa as visa
import Gpib
import math

class Agilent_MSO6034A: 
	
	def __init__(self, gpib): 
		
		###############################
		# Set Address and Reset Scope #
		###############################
		self.scope =  Gpib.Gpib(0, gpib)	
		self.scope.write("*IDN?")
		print self.scope.read(1024) 
		
		#############################
		# Initialize Scope Channels #
		#############################
		self.scope.write("CHAN1:COUP DC")
		self.scope.write("CHAN1:DISP 1")
		self.scope.write("CHAN2:COUP DC")
		self.scope.write("CHAN2:DISP 1")
		

		############################
		# Set up CH1:VOLT/CH2:AMPS #
		############################
		self.scope.write("CHAN1:UNIT VOLT")
		self.scope.write("CHAN2:UNIT AMP")
		#self.scope.write("CHAN2:OFFS 0")

		#############################
		# Set Vertical Scale Params #
		#############################
		self.scope.write("CHAN1:SCAL 2")
		self.scope.write("CHAN2:SCAL 0.020")

		###########################
		# Set up External Trigger #
		###########################
		self.scope.write("TRIG:SOUR EXT")
		self.scope.write("TRIG:SLOP POS")

		#############################
		# Set up Acquisition Params #
		#############################
		self.scope.write("ACQ:TYPE AVER")
		self.scope.write("ACQ:COUN 512")
		self.scope.write("STOP")

	def setTime(self, time):

		##################
		# Set Time Scale #
		##################
		self.scope.write("TIM:SCAL %sus"%(time))
		self.time = time

	def runMeas(self):
		self.scope.write("RUN")

	def stopMeas(self):
		self.scope.write("STOP")

	def readData(self):

		_dat = OrD()
		## Set the total number of points to be read
		self.scope.write("WAV:POIN 1000")

		## Read Channel 1
		self.scope.write("WAV:SOUR CHAN1")
		self.scope.write("WAV:FORM ASCII")
		self.scope.write("WAV:DATA?")
		_v = self.scope.read(16000)
		time.sleep(1)

		## Read Channel 2
		self.scope.write("WAV:SOUR CHAN2")
		self.scope.write("WAV:FORM ASCII")
		self.scope.write("WAV:DATA?")
		_i = self.scope.read(16000)
		time.sleep(1)
		
		## Get data into dictionary
		_dat["Voltage"] = [float(v) for v in re.split('[,]',_v)[1:-1]]
		_dat["Current"] = [float(i) for i in re.split('[,]',_i)[1:-1]]
	  
		## Voltage/Current Leveling and calcualte time axis
		_dat["Voltage"] = [v - _dat["Voltage"][0] for v in _dat["Voltage"]]
		_dat["Current"] = [i - _dat["Current"][0] for i in _dat["Current"]]
		_dat["Time"] = [self.time*n for n in range(len(_dat["Voltage"]))]
		return _dat