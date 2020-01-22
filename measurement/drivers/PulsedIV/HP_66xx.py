#!/usr/bin/env python 
#import pyvisa as visa
import Gpib

# Driver Code Here
class source_66xx:

	def __init__(self, GPIB):

		# Pyvisa
		# self.source = visa.instrument("GPIB::%s"%GPIB)
		
		# Linux GPIB
		self.source = Gpib.Gpib(0,gpib)

		# Query IDN
		self.source.ask("*IDN?")

		# Zero Voltages
		self.source.ask("VSET 1,0") 
		self.source.ask("VSET 2,0") 

		# Set outputs off
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
				self.source.write("OUT 1,0")
				self.source.ask("VSET 1,"+str(v1set)) 
				self.source.ask("ISET 1,"+str(compliance)) 
				self.source.write("OUT 1,1")
		
				# Output Two
				self.source.write("OUT 2,0")
				self.source.ask("VSET 2,"+str(v2set))
				self.source.ask("ISET 2,"+str(compliance)) 
				self.source.write("OUT 2,1")
		
			except: 
				print "Error: Check GPIB connection and address"
	
	def outputsOff(self):
	
		self.source.write("OUT 1,0")
		self.source.write("OUT 2,0")
	
	def getVoltage(self):  
	
		print "V1="+self.source.ask("VOUT? 1")
		print "V2="+self.source.ask("VOUT? 2")
	
	def getCurrent(self):
	
		print "I1="+self.source.ask("IOUT? 1")
		print "I2="+self.source.ask("IOUT? 2")