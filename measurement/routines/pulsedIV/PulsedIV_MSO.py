#!/usr/bin/env python 
import time

# Numpy and matplotlib
import numpy as np
import matplotlib.pyplot as plt

# Import insturments 
import AgilentMSO6034A
import source_66xx
import tabor_8500

# Measurement routine
if __name__ == "__main__": 

	## Vmax in VoltA and wPulse in us
	vmax, wPulse, measID = 8, 100, 'm4'

	## A generator function which provides 
	## the pulse voltages
	def genV(V):
		v = 1
		while v <= 1+4*V:
			_V = [n*0.25 for n in range(v)]
			yield _V+_V[::-1][1:]
			v += 4
   
	## Set up the generator to sweep to MAXIMUM 
	## voltage. It will generate the hysteresis 
	## of pulses up to vmax. Also skip the null 
	## list. 
	_v = genV(vmax)
	_v.next()
  
	## Get the Instruments on the go
	_scopegpib = 7
	scope = AgilentMSO6034A(_scopegpib)
	time.sleep(1)

	# Source 
	_sourcegpib = 2
	source = source_66xx(_sourcegpib)
	time.sleep(1)
   
	# Pulse 
	_pulsegpib = 9
	time.sleep(1)
	pulse = tabor_8500(_pulsegpib)
	pulse.setPulseA(wPulse,1)
	scope.setTime(wPulse)

	while True: 
		
		try:
			## Get the next vlist
			vlist = _v.next()

			## Set filename for data
			fname = './data/%s/%sus/%sV.pkl'%(measID,wPulse, int(max(vlist)))
		   
			## Run the measurement
			d = []
			for v in vlist: 

				print "Pulsing Voltage:%sV"%(v) 
				pulse.setVoltageA(0.0,5.0)
				time.sleep(2) 
				# Set Supply Voltage
				source.setVoltage(0,v,0.100)
 
				# Trigger Pulse/Scope
				pulse.triggerA()
		
				# Trigger Scope
				scope.runMeas()
				time.sleep(2) 
				scope.stopMeas()
				pulse.offPulseA()

				d.append(scope.readData())
				time.sleep(1) 
 
			# Dump the data
			
			print '--------------------------------------'
			print 'Dumping Data to %s'%(fname)
			print '--------------------------------------'
			pickle.dump(d, open(fname, "wb" ) )

		except StopIteration:
			break
