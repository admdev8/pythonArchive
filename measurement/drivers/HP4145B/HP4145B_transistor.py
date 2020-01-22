#!/usr/bin/env python 
import math
import visa
import re
import os 
import time

# Numerics and plotting
import ndfit
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mpl as mpl
import matplotlib.cm as cm

# HP4145B Semiconductor Parameter Analyzer
class hp_4145B:

	def __init__(self, GPIB): 
		## Initialize the dude
		self.analyzer = visa.instrument("GPIB::%s"%GPIB)

	def channelDefinition(self):
 
		## Select channel [DE]finition page
		self.analyzer.write("DE")
	   
		## Perform Channel Definition Steps. We need 
		## Three channels for a basic DC measurement.
		## 1 = Voltage Mode
		## 2 = Current Mode
		## 3 = Common
		## last value = 1 ---> VAR1 is active

		self.analyzer.write("DE CH1,'V1','I1',1,1") ## Drain/Gate
		self.analyzer.write("DE CH4,'V4','I4',1,2") ## Gate/Drain
		
		self.analyzer.write("DE CH2,'VB','IB',3,3") ## Back Gate Common
		self.analyzer.write("DE CH3,'VS','IS',3,3") ## Source Common
		
		## Turn off other channels 
		self.analyzer.write("DE VS1")
		self.analyzer.write("DE VS2")
		self.analyzer.write("DE VM1")
		self.analyzer.write("DE VM2")


	def sweepVoltage(self,start, stop, step,compliance):
		## Select [S]weep [S]etup page  
		self.analyzer.write("SS")

		## Set up the gate voltage variable. 
		CMD = "VR1,%s,%s,%s,%s"%(str(start), str(stop), str(step), str(compliance))
		self.analyzer.write(CMD)
		self.number = np.abs(round((stop-start)/step,1))+1

		## Return a list of the sweep variables
		return list(np.linspace(start,stop,self.number))
	   
	def stepVoltage(self,start,step,number,compliance):
			## Set up drain voltage setup
		self.analyzer.write("SS")
		CMD = "VP %s,%s,%s,%s"%(str(start), str(step), str(number), str(compliance))
		self.analyzer.write(CMD)

	def setOutput(self):
		## Set up the list display. Note that if something is not in
		## the list we will not be able to request the data.
		self.analyzer.write("SM")
		self.analyzer.write("DM2")
		self.analyzer.write("LI 'I1','I4'")

		## Pack the independent parameter into a list and return it
		#current = list(np.linspace(start, stop, int((stop-start)/(step)), endpoint=False))
		#current.append(stop)
		#return current

	def measureData(self):
		## Issue command to [M]easure [D]ata
		self.analyzer.write("MD ME1")

		## Pause for the 4145B to populate the buffer with measurement
		## data. If there is no sleep, Python will issue the ask commend
		## before the measurement has been complete

	def getDataIdVd(self):
		time.sleep(5)
		## Request data and transform it in to a list by splitting 
		## on regular expression. Note that the slice operation is 
		## needed because the first item in the list after re.split
		## is the empty string.
		data = re.split(',*[NXCT]\s*', self.analyzer.ask("DO 'I1'"))[1:] 
		return [float(i) for i in data]
	
	def getDataIdVg(self):
		time.sleep(5)
		## Request data and transform it in to a list by splitting 
		## on regular expression. Note that the slice operation is 
		## needed because the first item in the list after re.split
		## is the empty string.
		data = re.split(',*[NXCT]\s*', self.analyzer.ask("DO 'I4'"))[1:] 
		return [float(i) for i in data]


## Useful Generator Function 
def chunks(l, n):
    for i in xrange(0, len(l), n):
        yield l[i:i+n]


## Mainloop
if __name__ == "__main__": 
	
	## Gate Width
	w = 10e-6
	
	## For Vds/Ids measurements
   
	#############################
	##		SMU1 DRAIN	   ##
	##		SMU4 GATE		##
	#############################
	
	## Cheap Switch
	if 1 == 1: 
		## Some Gate Voltages
		vgx = list(np.linspace(0,0.8,9)) 
		## A program to test things out
		GPIB = 1
		pa = hp_4145B(GPIB)
		pa.channelDefinition()
		
		## Sweep Vd from 0V to 1.1V
		Vd = pa.sweepVoltage(0.0,1.0,0.05,"40E-3")

		fig = plt.figure(figsize=(10,8))
		ax1 = fig.add_axes([0.1, 0.1, 0.70, 0.8])
		ax2 = ax1.twinx()
		for n,vg in enumerate(vgx):
			# time.sleep(1)
			pa.stepVoltage(vg,0.0,1,"30E-3")
			pa.setOutput()
			pa.measureData()
			Id = pa.getDataIdVd()
			g0 = [xx/w for xx in  list(ndfit.derivative(Id,Vd))]
		 
			print "%s %f"%("Vg =", vg)

			c = cm.cool(float(n)/len(vgx),1)
			ax1.plot(Vd, Id,color=c)
			ax2.plot(Vd, g0,color=c,linestyle='--')

		ax1.set_ylim(0.0, max(Id)+0.001)
		ax2.set_ylim(0.0, max(g0)+100)
		ax1.set_xlabel(r"Drain Voltage $(V_d)$")
		ax1.set_ylabel(r"Drain Current $(I_d)$")
		ax2.set_ylabel(r"Output Conductance $(g_0)$")

		## Fancy Colorbar
		cmap = mpl.cm.cool
		norm = mpl.colors.Normalize(min(vgx),max(vgx))
		ax2 = fig.add_axes([0.90, 0.1, 0.03, 0.8])
		cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,
										norm=norm,
										orientation='vertical')

		cb1.set_label(r"Gate Voltage $(V_g)$")
		plt.show()

	## For Vgs/Ids measurements
   
	#############################
	##		SMU1 GATE		##
	##		SMU4 DRAIN	   ##
	#############################

	if 1 == 0: 
		
		## Some Drain Voltages
		vdx = list(np.linspace(0,0.8,8)) 
		## A program to test things out
		GPIB = 1
		pa = hp_4145B(GPIB)
		pa.channelDefinition()
		
		## Sweep Vg from 0V to 1V
		Vg = pa.sweepVoltage(0.0,0.9,0.05,"1E-6")
		
		fig = plt.figure(figsize=(10,8))
		ax1 = fig.add_axes([0.1, 0.1, 0.70, 0.8])
		ax2 = ax1.twinx()
		for n,vd in enumerate(vdx):
			# time.sleep(1)
			pa.stepVoltage(vd,0.0,1,"20E-3")
			pa.setOutput()
			pa.measureData()
			Id = pa.getDataIdVg()
   
			print "%s %f"%("Vd =", vd)
			gm = [xx/w for xx in  list(ndfit.derivative(Id,Vg))]

			c = cm.cool(float(n)/len(vdx),1)
			ax1.plot(Vg, Id,color=c)
			ax2.plot(Vg,gm,color=c,linestyle='--')


		ax1.set_ylim(0.0, max(Id)+0.001)
		ax2.set_ylim(0.0, max(gm)+50)
		ax1.set_xlabel(r"Gate Voltage $(V_g)$")
		ax1.set_ylabel(r"Drain Current $(I_d)$")
		ax2.set_ylabel(r"Transconductance $(g_m)$")
		
		
		## Fancy Colorbar
		cmap = mpl.cm.cool
		norm = mpl.colors.Normalize(min(vdx),max(vdx))
		ax2 = fig.add_axes([0.90, 0.1, 0.03, 0.8])
		cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,
										norm=norm,
										orientation='vertical')

		cb1.set_label('Drain Voltage $(V_d)$')
		plt.show()
