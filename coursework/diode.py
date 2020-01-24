#!/usr/bin/env python 
from matplotlib import pyplot as plt
import numpy as np
import ndfit
import sys

## For Plotting
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':20})
rc('text', usetex=True)
rc('xtick', labelsize=20) 
rc('ytick', labelsize=20) 

## Set max recursion depth 
sys.setrecursionlimit(1500) 

## Constants of the Problem (global)
G,l,Is,sx = 0.3333,29.585,1e-11,3
class diode:

	## Some values to store the number of steps
	def __init__(self):
		self.convergence = [] 
		self.counter = 0

	## A few equations for the diode model
	def _diode(self,v): return Is*(np.exp(l*v)-1)
	def ddiode(self,v): return l*Is*(np.exp(l*v))

	## High bias model from [A]
	def xdiode(self,v): return Is*(np.exp(l*sx)*(1 + l*(v-sx))) - Is
	def xddiode(self,v): return l*Is*(np.exp(l*sx))

	def check(self,v,Es): return self._diode(v) - G*Es + G*v

	## Newtons method (tail recursive - value storing)
	def solve(self,v, Es):

		## Set a counter
		self.counter+=1
		
		## Add to the list to monitor convergence
		self.convergence.append(v)

		## Do the calculation  
		if (v < sx):
			vm = v - (self._diode(v)- G*(Es-v))/(self.ddiode(v)+G)
		else:
			vm = v - (self.xdiode(v)- G*(Es-v))/(self.xddiode(v)+G)


		## Check convergence 0.01mV
		if self.check(vm, Es) < .001:
			self.counter,self.check(vm,Es)
			self.convergence.append(vm)
			self.result = vm
			return 

		elif self.counter > 1500:
			self.result = vm
			return 
		   
		## Make recursive call
		else:
			self.counter+=1
			self.convergence.append(vm)
			self.solve(vm,Es)
		
## main loop 
if __name__ == "__main__":

	## Generate a list of voltages
	_Es = np.linspace(-.5,3,100)
	_Si = [3*np.sin(2*i) for i in np.linspace(0,2*3.14,200)]
   
	## Solve Newton for everybody in the list
	result,sresult = [],[]
	for i in _Es:
		d = diode()
		d.solve(.5*i,i)
		result.append(d.result)

	## And for the Sine wave
	for i in _Si:
		d = diode()
		d.solve(i*.5,i)
		sresult.append(d.result)

	## Plot of the operating point vs. bias point  
	fig = plt.figure(1)
	ax1 = fig.add_subplot(111)
	ax2 = ax1.twinx()
	ax1.plot(_Es,result)
	ax2.plot(_Es,ndfit.quotient(list(result), list(_Es)),'r')
	ax1.set_xlabel(r"Applied Bias $(E_s)$")
	ax1.set_ylabel(r"Diode Voltage $(V_D)$")
	ax2.set_ylabel(r"Voltage Ratio $(V_D/E_s)$")
	ax1.set_title(r"Diode Operating Point")

	## Plot of waveform response
	fig = plt.figure(2)
	ax1 = fig.add_subplot(111)
	ax2 = ax1.twinx()
	ax1.plot(sresult)
	ax1.set_ylim(-3,3)
	ax2.plot(_Si,'--r')
	ax1.set_xlabel(r"Time Step")
	ax1.set_ylabel(r"Diode Voltage $(V_D)$")
	ax2.set_ylabel(r"Input Waveform $(E_s)$")
	ax1.set_title(r"Waveform Response")
	
	## Plot of convergence at Es = 1.2V
	d = diode()
	d.solve(0.6,1.2)
	cur = [d._diode(i) for i in d.convergence]
	fig = plt.figure(3)
	ax1 = fig.add_subplot(111)
	ax2 = ax1.twinx()
	ax1.plot(d.convergence)
	ax2.semilogy(cur,'r')
	ax1.set_xlabel(r"Iteration")
	ax1.set_ylabel(r"Diode Voltage $(V_D)$")
	ax2.set_ylabel(r"Diode Current $(I_D)$")
	ax1.set_title(r"Diode Convergence (1.2V)")

	## Plot of convergence at Es = 15V
	d = diode()
	d.solve(15,15)
	cur = [d._diode(i) for i in d.convergence]
	fig = plt.figure(4)
	ax1 = fig.add_subplot(111)
	ax2 = ax1.twinx()
	ax1.plot(d.convergence)
	ax2.semilogy(cur,'r')
	ax1.set_xlabel(r"Iteration")
	ax1.set_ylabel(r"Diode Voltage $(V_D)$")
	ax2.set_ylabel(r"Diode Current $(I_D)$")
	ax1.set_title(r"Diode Convergence (15V)")
	plt.show()
