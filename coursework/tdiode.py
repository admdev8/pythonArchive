from matplotlib import pyplot as plt
import matplotlib.mpl as mpl
import matplotlib.cm as cm
import numpy as np
import ndfit
import sys

## For timing siulations 
import timeit

## For Plotting
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':20})
rc('text', usetex=True)
rc('xtick', labelsize=20) 
rc('ytick', labelsize=20)
rc('legend', fontsize=20) 

## Set max recursion depth 
sys.setrecursionlimit(15000)

## Make a diode class. This stores our diode model
class diode:

	## Allow the resistnace to be changed
	def __init__(self):
		## Constants of the problem
		self.vt,self.n,self.Is,self.Cj,self.phi,self.g = 0.026,1.3,1e-11,1.2e-12,1.3,0.5

	## The diode capacitance
	def diodec(self,vd): return self.Cj * np.power((1-(vd/self.phi)),-self.g)

	## The diode current
	def diodei(self,vd): return self.Is*(np.exp(vd/self.n/self.vt)-1)


class tdiode:

	## Allow the resistnace to be changed
	def __init__(self):
		## Constants of the problem
		self.vt,self.n,self.Is,self.Cj,self.phi,self.g = 0.026,1.3,1e-11,1.2e-12,1.3,0.5
		self.r0,self.m,self.v0 = 1, 2.0, 0.3

	## The diode current
	def diodei(self,vd):

		## The tunnel current
		It = (vd/self.r0)*np.exp(-np.power((vd/self.v0),self.m))

		## Ordinary diode current
		Id = self.Is*(np.exp(vd/self.n/self.vt)-1)

		return It + Id 
	
	## The diode capacitance
	def diodec(self,vd): return self.Cj * np.power((1-(vd/self.phi)),-self.g)


## A class to calculate the transient of
## the diode for various series resistances
class transient:

	## Class initializer
	def __init__(self, V, delta, G):

		## Save the input waveform and series resistance 
		self.V, self.G,self.delta = V,G,delta 
		
		## Create a diode object by calling diode constructor
		self.d = tdiode()
		self.result, self.counter = [],0

		## Lambdas for backwards iteration and trapeziodal rule
		self.b = lambda vn,vm,V: vn - vm - (self.circuit(V,vn)*self.delta)
		self.t = lambda vn,vm,V: vn - vm - (0.5*self.delta*(self.circuit(V,vn) + (self.circuit(V,vm))))

	## A method to calculate the step for our circuit
	def circuit(self,V, vm):
		return (self.G*(V-vm) - self.d.diodei(vm))/self.d.diodec(vm)

	## Forward Euler calculation (recursive)
	def forward(self, vm):
		try:
			## Calculate next step
			vn = vm + self.circuit(self.V[self.counter],vm)*self.delta
			self.result.append(vn)
			self.counter+=1
  
			## Make recursive call
			self.forward(vn)
			
		except IndexError:
			return  

	## Forward Euler calculation (recursive)
	def backward(self, vm):
		try:
			## Calculate next step
			vn = _solve(self.b, vm, self.V[self.counter])
			self.result.append(vn)
			self.counter+=1
  
			## Make recursive call
			self.backward(vn)
			
		except IndexError:
			return  

	## Forward Euler calculation (recursive)
	def trapezoid(self, vm):
		try:
			## Calculate next step
			vn = _solve(self.t, vm, self.V[self.counter])
			self.result.append(vn)
			self.counter+=1
  
			## Make recursive call
			self.trapezoid(vn)
			
		except IndexError:
			return 

## We need to find vn given vm (Newton Raphson/iterative)
conv,d = 1e-6,1e-6
def _solve(func,_vm,_V,_vn=0.0):

	while True:
		r = func(_vn,_vm,_V)
		if np.abs(r) < conv:
			return _vn
		else:
			dfunc = (func(_vn+d,_vm,_V) - r)/d
			_vn = (_vn - (r/dfunc))
			r = func(_vn,_vm,_V)


## Make a squarepulse
def squarepulse(low, high, time, step):
	npoints = int(time/step)
	wave = [low]*npoints + [high]*npoints + [low]*npoints + [low]*npoints 
	time = [i*step for i in range(4*npoints)]
	return time, wave, step
	
if __name__ == "__main__":


	if 1 == 0:

		## Initialize the diode
		t = tdiode()

		## A list of operating points
		vd = np.linspace(0,0.80, 100)

		## Calculate the current
		i = [t.diodei(_vd) for _vd in vd]

		plt.plot(vd,i)
		plt.show()

	if 1 == 0: 

		## Generate initial waveform for plotting
		time,wave,step = squarepulse(0.0,0.4,3e-9,1e-12)

		## Plot the results
		fig = plt.figure(figsize=(10,8))
		ax1 = fig.add_axes([0.1, 0.1, 0.7, 0.8])
		ax1.plot(time, wave, "k--", label="$V(t)$")
		l = ax1.legend()
		l.draw_frame(True) 
		ax1.set_xlabel('Time $(s)$')
		ax1.set_ylabel('Diode Voltage $(V_D)$')

		## Set up a list of time steps
		tlist = list(np.linspace(1e-12,4e-10,8))
		for i,t in enumerate(tlist[::-1]):
			time,wave,step = squarepulse(0.0,1.0,3e-9,t)
			c = cm.cool(float(i)/len(tlist),1)
			f = transient(wave,step,0.005)
			f.trapezoid(0.0)
			ax1.plot(time,f.result,color= c)
			ax1.set_ylim(0,1.1)

		ax1.set_title("Trapezoidal")
		## Fancy colorbar 
		cmap = mpl.cm.cool_r
		norm = mpl.colors.Normalize(min(tlist)/1e-12,max(tlist)/1e-12)
		ax2 = fig.add_axes([0.85, 0.1, 0.03, 0.8])
		cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, orientation='vertical')
		cb1.set_label('Time Step $(ps)$')
		plt.show()

	## A cheap switch to turn on/off this simulation 
	if 1 == 1: 

		
		## Initialize the diode
		t = tdiode()

		## A list of operating points
		vd = np.linspace(0,0.80, 100)

		## Calculate the current
		i = [t.diodei(_vd) for _vd in vd]

		plt.plot(vd,i)

		######################################
		######################################
		
		tic = timeit.default_timer()
		
		## Initialize empty list for simulation results
		simlist = []
		time,wave,step = squarepulse(0.0,1.0,3e-9,1e-12)
		
		## Generate a list of input impedances
		rlist = np.linspace(1,10,20)
		glist = [1./r for r in rlist]

		## Simulate transient for each input impedance
		for g in glist:
			t = transient(wave,step,g) 
			t.backward(0.0)
			simlist.append(t.result)

		toc = timeit.default_timer()
		print toc - tic 
			
		## Plot the results
		fig = plt.figure(figsize=(10,8))
		ax1 = fig.add_axes([0.1, 0.1, 0.7, 0.8])
		ax1.plot(time, wave, "k--", label="$V(t)$")
		l = ax1.legend()
		l.draw_frame(True) 
		ax1.set_xlabel('Time $(s)$')
		ax1.set_ylabel('Diode Voltage $(V_D)$')
			
		for i,result in enumerate(simlist):
			c = cm.cool(float(i)/len(rlist),1)
			ax1.plot(time,result,color=c)
			ax1.set_ylim(0,1.1)

		## Fancy colorbar 
		cmap = mpl.cm.cool
		norm = mpl.colors.Normalize(min(rlist),max(rlist))
		ax2 = fig.add_axes([0.85, 0.1, 0.03, 0.8])
		cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, orientation='vertical')
		cb1.set_label('Source Impedance $(\Omega)$')
		plt.show()
