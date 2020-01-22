#!/usr/bin/env python 
#import pyvisa as visa
import Gpib
import math

class Tektronix_TDS7104:

	def __init__(self, GPIB):

		self.scope = visa.instrument("GPIB::%s"%GPIB)

		#############################
		# Initialize Scope Channels #
		#############################
		
		# Set Horizontal Axis Parameters
		self.scope.write("HOR:SCA 1E-6")
		self.scope.write("HOR:TRIG:POS: 20%")
		# Set Horizontal Resolution 
		self.scope.write("HOR:RESO 5000")
		self.scope.write("HOR:RECO 5000")

		# Zero Vertical Offsets of Channels to Zero
		self.scope.write("CH1:OFF 0")
		self.scope.write("CH2:OFF 0")
		# Set channel copuling to DC 
		self.scope.write("CH1:COUP DC")
		self.scope.write("CH2:COUP DC")
		# Set Reasonable Vertical Scale of Channels
		self.scope.write("CH1:SCA 2")
		self.scope.write("CH2:SCA 0.02")
		self.scope.write("CH3:SCA 2")
		self.scope.write("CH4:SCA 0.02")
		
		# Turn off acquisition
		self.scope.write("ACQ:STATE OFF")
		self.scope.write("ACQ:REPE OFF")
				
		#####################
		# Set up Triggering #
		#####################

		# Trigger from AUX in. Rising Slope with DC coupling
		self.scope.write("TRIG:A:EDGE:COUP DC")
		self.scope.write("TRIG:A:EDGE:SLO RIS")
		self.scope.write("TRIG:A:EDGE:SOU AUX")

	def setScale(self, channel, scale): 
		self.scope.write("CH%s:SCA %s"%(channel, scale)) 
	
	def setSingle(self):
		# Enable Single Aquisition Mode
		self.scope.write("ACQ:STATE ON")
		self.scope.write("ACQ:REPE ON")
		self.scope.write("ACQ:STOPA SEQ")

	def setAverage(self, numavg): 
		# Enable averaging mode
		self.scope.write("ACQ:MOD AVE")
		self.scope.write("ACQ:STOPA SEQ")
		self.scope.write("ACQ:NUMAV %s"%numavg)
		self.scope.write("ACQ:STATE ON")

	def setAcqOff(self): 
		# Write to turn off aquistion this prevents 
		# averaging of sucsessive aquisitions
		self.scope.write("ACQ:STATE OFF")

	def chx_data(self, channel):
		
		# Request Encoding and curve Data
		self.scope.write("DAT:INIT")
		self.scope.write("DAT:SOU CH%s"%channel)
		self.scope.write("DAT:ENC ASCIi")
		datString = self.scope.ask("CURV?")
		
		# Obtain x and y scaling factors
		nTstep = int(self.scope.ask("HOR:RECO?"))
		tScale = float(self.scope.ask("WFMO:XIN?"))
		vScale = float(self.scope.ask("WFMO:YMU?"))

		# Build Time Data
		tList = range(0,nTstep) 
		timeData = [tScale*float(i) for i in tList] 

		# Build Votlage Data
		my_splitter = shlex.shlex(datString, posix=True)
		my_splitter.whitespace += ','
		my_splitter.whitespace_split = True
		vList = list(my_splitter)
		voltageData = [vScale*float(i) for i in vList]

		return timeData, [i-voltageData[0] for i in voltageData]
		