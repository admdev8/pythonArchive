#!/usr/bin/env python 
#import pyvisa as visa
import Gpib
import math

# Class for Tabor 8500 pulse generator
class tabor_8500: 
 
	def __init__(self, gpib):

		# Set Address 
		# self.pulse = visa.instrument("GPIB::%s"%GPIB)

		# Set Address 
		self.pulse = Gpib.Gpib(0, gpib)	
		
		# Select Channel A and turn output off
		self.pulse.write("CHA")
		self.pulse.write("D1") 
				
	def setPulseA(self,wPulse,Period):  
		
		self.time = wPulse
		
		# Set Parameters
		self.pulse.write("WID %sus"%wPulse) 
		self.pulse.write("PER %sms"%Period) 

		# Activate Output
		self.pulse.write("D0")
	
	def setVoltageA(self,vLow,vHigh):

		self.vLow, self.vHigh = vLow, vHigh

		if math.fabs(vHigh) > 5:

			print "Voltage is too high: output off"
			self.vHigh = 0
		
		if math.fabs(vLow) > 5:

			print "Voltage is too high: output off"
			self.vLow = 0
	
	def offPulseA(self):
	
		self.pulse.write("D1")
		
	def triggerA(self):
		
		if self.vLow != self.vHigh: 
			
			self.pulse.write("LOL %sV"%round(self.vLow,2))
			self.pulse.write("HIL %sV"%round(self.vHigh,2))
			
			# Set to Triggered Mode and turn on output
			self.pulse.write("D0")
			self.pulse.write("M1")
			self.pulse.write("TRG")

		else: 
			self.pulse.write("D1")
