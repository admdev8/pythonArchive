#!/usr/bin/python2.6
import math
import visa
import shlex
import time
import os

import matplotlib.pyplot as plt
# Connect to Devices

#########
# Usage #
#########

if __name__ == "__main__":

	# Power Supply
	GPIB, vlow, vhigh, compliance = 6,0,4,.5
	source = source_66xx(GPIB)
	source.setVoltage(vlow,vhigh,compliance) 
	source.getVoltage() 
	source.getCurrent()

	# Scope
	GPIB = 15
	scope = scope_tds7104(GPIB)
	scope.setAverage(512)
	scope.setScale(1,2) 
	scope.setScale(2,0.02) 
	scope.setScale(3,2)

	# Pulse 
	GPIB = 9 
	pulse = tabor_8500(GPIB)
	pulse.setPulseA(1,1)
	pulse.setVoltageA(0,2)
	pulse.triggerA()

	# Buffer sleep
	time.sleep(1)

	# Scope Data
	tch1, vch1 = scope.chx_data(1)
	tch2, vch2 = scope.chx_data(2)
	tch3, vch3 = scope.chx_data(3)

	plt.plot(tch1,vch1)
	plt.plot(tch2,vch2)
	plt.plot(tch3,vch3)
	plt.title("Test Pulses")
	plt.xlabel("Time (s)")
	plt.ylabel("Voltage (V)")
	plt.show()
	#source.outputsOff()
