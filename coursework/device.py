#!/usr/bin/env python 
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.mpl as mpl
import numpy.linalg as la
import matplotlib.cm as cm
from numpy import unravel_index
import numpy as np
import copy
import math
import sys
import os

## Custom modules
import microwaveConverter as mwc
import minispice as ms
import ndfit

## For timing simulations 
import timeit


## For Plotting
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':20})
rc('text', usetex=True)
rc('xtick', labelsize=20) 
rc('ytick', labelsize=20)
rc('legend', fontsize=20) 

## Try 3D plotting
import pylab as p
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3

## Set max recursion depth 
sys.setrecursionlimit(15000)

## List of constants for model
a, b, c, d, g, phi = -0.080, 3.3, 0.41, 10.0, 0.056, math.pi*2.4/180

def _uplus(Vgd,Vgs): 
	pgd =  np.cos(phi)*Vgd + np.sin(phi)*Vgs 
	pgs =  np.cos(phi)*Vgs - np.sin(phi)*Vgd 
	return pgd, pgs

def _uminus(Vgd,Vgs): 
	mgd =  np.cos(phi)*Vgd - np.sin(phi)*Vgs 
	mgs =  np.sin(phi)*Vgd + np.cos(phi)*Vgs 
	return mgd, mgs

## NEED f1(ugd+ , ugs+)*f2(vgs-vgd) - f1(ugs- , ugd-)*f2(vgd-vgs)
def f1(u1,u2):
	const = np.exp(-b*(u2+c))
	return (1+a*u1)*(1 - np.tanh(const))

def f2(v):
	const = np.exp(-d*v)
	return (1 - np.tanh(const))
		
def current(Vgd, Vgs):
	pgd, pgs = _uplus(Vgd, Vgs) 
	mgd, mgs = _uminus(Vgd, Vgs) 
	delta = Vgs - Vgd
	return  g*(f1(pgd, pgs)*f2(delta) - f1(mgs,mgd)*f2(-delta))


##################################
## EXTREMELY VALUABLE FUNCTIONS ##
##################################
def plotcontour(x,y,z,xlabel=None, ylabel=None, title=None, arg=None):
   ## Plot the first deriv Ix
   fig = plt.figure(figsize=(10,8))
   ax1 = fig.add_axes([0.1, 0.1, 0.7, 0.8])

   ## Fancy colorbar 
   cmap = mpl.cm.RdBu
   norm = mpl.colors.Normalize(np.amin(z),np.amax(z))

   if arg is not None: 
	   ax1.contourf(x,y,z,99,cmap=cm.RdBu)
	   ax1.set_xscale('log')
	   ax1.set_yscale('log')
   else:
	   ax1.contourf(x,y,z,30,cmap=cm.RdBu)

   ax1.set_xlabel(xlabel)
   ax1.set_ylabel(ylabel)
   ax2 = fig.add_axes([0.85, 0.1, 0.03, 0.8])
   cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,
								   norm=norm,
								   orientation='vertical')
   cb1.set_label(title)

#####################################
#		DERIVATIVE FUNCTIONS	   #
#####################################
def _partialx(Z,X):
	return np.array([ndfit.derivative(list(Z[i,:]), list(X[i,:])) for i in range(len(Z))])
def _partialy(Z,Y):
	return np.array([ndfit.derivative(list(Z[:,i]), list(Y[:,i])) for i in range(len(Z))]).T

def geti1(h,pd,a):
	return (a)*getg1(h,pd)
def geti2(h,pd,a):
	return (a**2)*getg2(h,pd)
def geti3(h,pd,a):
	return (a**3)*getg3(h,pd)
def getg1(h,pd):
	return (pd['d']*(1-h) + pd['s'])
def getg2(h,pd):
	return (1./2)*(pd['dd']*((1-h)**2) + 2*pd['ds']*(1-h) + pd['ss'])
def getg3(h,pd):
	return (1./6)*(pd['ddd']*((1-h)**3) + 3*pd['dds']*((1-h)**2) + 3*pd['dss']*(1-h) + pd['sss'])


#####################################
#   SPECIFIC SIMULATION EXAMPLES	#
#####################################
if __name__ == "__main__":

	if 1 == 0:
		vgs = np.arange(-3, 0.5, 0.01)
		vgd = np.arange(-3, 0.5, 0.01)
		_vgd,_vgs = np.meshgrid(vgd, vgs)
   
		## Get index of bias point (dumb way)
		VGD,VGS = -2.55,-0.55  
		for row,val in enumerate(vgd):
			if round(val,2)==VGD:
				print val
				break

		for col,val in enumerate(vgs):
			if round(val,2)==VGS:
				print val
				break

		######################################
		## CALCULATE CURRENT AND ITS DERIVS ##
		######################################

		## Get the current
		Ids = current(_vgd,_vgs)

		## Get first derivatives 
		Ix = _partialx(Ids,_vgd)
		Iy = _partialy(Ids,_vgs)

		## Get the second derivatives
		Ixx = _partialx(Ix,_vgd)
		Ixy = _partialy(Ix,_vgs)
		Iyy = _partialy(Iy,_vgs)

		## Get the third Derivatives
		Ixxx = _partialx(Ixx,_vgd)
		Ixxy = _partialy(Ixx,_vgs)
		Iyyx = _partialy(Ixy,_vgs)
		Iyyy = _partialy(Iyy,_vgs)

		## Create a dictionary of the partials
		## at the bias point.
		pd = {}
   
		## The constant term
		pd[''] =Ids[col,row]

		## The first derivs
		pd['d']  = Ix[col,row]
		pd['s']  = Iy[col,row]

		## The second derivs
		pd['dd'] = Ixx[col,row]
		pd['ds'] = Ixy[col,row]
		pd['ss'] = Iyy[col,row]

		## The third derivs
		pd['ddd']= Ixxx[col,row]
		pd['dds']= Ixxy[col,row]
		pd['dss']= Iyyx[col,row]
		pd['sss']= Iyyy[col,row]

	################################
	#	  Linear Amplifier		#
	################################
	if 1==1:
		freqList = np.logspace(9.5, 10.5, 100)
		path=os.getcwd()+"/amplifier.cir"
		data = ms.frequencyAnalysis.fromFile(path, freqList)
		data.plotGain(1,9,"logdb")
		data.plotTransferFunction(1,9,"log")
		data.inOutImpedance(2,8,50,50,"log")
		plt.show()

	if 1==0:
		freqList = [1e10]
		path=os.getcwd()+"/amplifier.cir"
		data = ms.frequencyAnalysis.fromFile(path, freqList)
		yint =data.ygroup[0].toTwoport(2,8)

	##########################################
	#  Fancy plots of I and its derivatives  # 
	##########################################
	if 1==0:
		##  The Current Contour
		plotcontour(vgd,vgs,Ids,'$I_{ds}$ $(A)$')

		## First Derivatives
		plotcontour(vgd,vgs,Ix,'$\partial_{d}I$')
		plotcontour(vgd,vgs,Iy,'$\partial_{s}I$')

		## Second Derivatives
		plotcontour(vgd,vgs,Ixx,'$\partial_{dd}I$')
		plotcontour(vgd,vgs,Ixy,'$\partial_{ds}I$')
		plotcontour(vgd,vgs,Iyy,'$\partial_{ss}I$')
 
		## Third Derivatives
		plotcontour(vgd,vgs,Ixxx,'$\partial_{ddd}I$')
		plotcontour(vgd,vgs,Ixxy,'$\partial_{dds}I$')
		plotcontour(vgd,vgs,Iyyx,'$\partial_{dss}I$')
		plotcontour(vgd,vgs,Iyyy,'$\partial_{sss}I$')
		plt.show()
	  
	####################################
	#		 PLOT IV CURVES		   #
	####################################
	if 1 == 0:
		## Set up vgd and vds
		vgs, vds = np.linspace(-1, 1, 100), np.linspace(0,3,100)

		Ivec = []
		for _vgs in vgs: 
			I = []
			for _vds in vds: 
				I.append(current((_vgs-_vds),_vgs))
			Ivec.append(I)

			
		## Plot the results
		fig = plt.figure(figsize=(10,8))
		ax1 = fig.add_axes([0.1, 0.1, 0.7, 0.8])
		ax1.set_xlabel('$V_{ds}$ $(V)$')
		ax1.set_ylabel('$I_{ds}$ $(mA)$')
	
		for i,result in enumerate(Ivec):
			c = cm.cool(float(i)/len(vgs),1)
			ax1.plot(vds,result,color=c)
		   
		## Fancy colorbar 
		cmap = mpl.cm.cool
		norm = mpl.colors.Normalize(min(vgs),max(vgs))
		ax2 = fig.add_axes([0.85, 0.1, 0.03, 0.8])
		cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,
										norm=norm,
										orientation='vertical')

		cb1.set_label('$V_{gs}$ $(V)$')
		plt.show()
		

	################################
	#	  IP2/IP3 vs Vgs/Vds	  #
	################################
	if 1==0:

		npoints = 50
		vgs = np.linspace(0.7,-0.7,npoints)
		vds = np.linspace(0.0,3.0,npoints)
		_vgs,_vds = np.meshgrid(vgs,vds)
		_vgd = np.subtract(_vgs,_vds)
	   
		## Get the current
		_ids = np.zeros([npoints,npoints])
		for i,VGS in enumerate(vgs):
			for j,VDS in enumerate(vds):
				VGD = VGS-VDS 
				_ids[j,i] = current(VGD,VGS)

		## Test our initial bias points
		Vds,Vgs = 2.0,-0.55
		for row,val in enumerate(vds):
			if round(val,2)==Vds:
				break
		for col,val in enumerate(vgs):
			if round(val,3)==Vgs:
				break

		print Vgs, Vds
			
		## Now we need the effective Gm at each bias point.
		## Generally we calculate the derivitives.
		Ix = _partialx(_ids,_vgs)
		Iy = _partialy(_ids,_vgd)

		## ## Get the second derivatives
		Ixx = _partialx(Ix,_vgs)
		Ixy = _partialx(Iy,_vgs)
		Iyy = _partialy(Iy,_vgd)
		
		## ## Get the third Derivatives
		Ixxx = _partialx(Ixx,_vgs)
		Ixxy = _partialx(Ixy,_vgs)
		Iyyx = _partialx(Iyy,_vgs)
		Iyyy = _partialy(Iyy,_vgd)

		## Create a dictionary of the partials
		## at the bias point.
		pd2 = {}
   
		## The constant term
		pd2[''] =_ids[row,col]

		## The first derivs
		pd2['s']  = Ix[row,col]
		pd2['d']  = Iy[row,col]
	   
		## The second derivs
		pd2['ss'] = Ixx[row,col]
		pd2['ds'] = Ixy[row,col]
		pd2['dd'] = Iyy[row,col]

		## The third derivs
		pd2['sss']= Ixxx[row,col]
		pd2['dss']= Ixxy[row,col]
		pd2['dds']= Iyyx[row,col]
		pd2['ddd']= Iyyy[row,col]
		
		# Get the Gm and Gds
		_Gm  = np.add(Ix,Iy)
		_Rds = np.zeros([npoints,npoints])
		for i in range(len(Iy)):
			for j in range(len(Iy)):
				_Rds[i,j] = -1./Iy[i,j]

		## Set up a basic frequency analysis and initialize
		## the elements to the initial values
		f1,f2,f3 = 1e10,1e10+1e5,1e10+2e5
		freqList = [f1,f2,f3]
		path=os.getcwd()+"/amplifier.cir"
		data = ms.frequencyAnalysis.fromFile(path, freqList)

		## Set Rs, Rl, Gm, Rds to those of the original amplifier
		## in case they have been modified.
		Rs,Rl,Gm,Rds = 50.0,50.0,0.080348,228.653  
		data.setElement('Rs',Rs)
		data.setElement('Rl',Rl)
		data.setElement('Gm',Gm)
		data.setElement('Rds',Rds)

		## Now Proceed with the calculation of IP3. First we need some
		## Arrays to store the IP2 and IP3 values
		_IP2,_IP3 = np.zeros([npoints,npoints]),np.zeros([npoints,npoints])
		_GH,_GP = np.zeros([npoints,npoints]),np.zeros([npoints,npoints])
		
		## Now loop through all of the gm and rds.
		## calculate IP2 and IP3
		_v0 = np.logspace(-4,2,300)
		for (i,j), value in np.ndenumerate(_Gm):

			## Write the values to the circuit file.
			## Perform a frequency analysis
			data.setElement('Gm',_Gm[i,j])
			data.setElement('Rds',_Rds[i,j])
			data = ms.frequencyAnalysis.fromFile(path, freqList)

			## Get the correct value out of all the partials
			## These have already been evaluated at all bias points
			pd2 = {}
   
			## The constant term
			pd2[''] =_ids[i,j]

			## The first derivs
			pd2['s']  = Ix[i,j]
			pd2['d']  = Iy[i,j]
	   
			## The second derivs
			pd2['ss'] = Ixx[i,j]
			pd2['ds'] = Ixy[i,j]
			pd2['dd'] = Iyy[i,j]

			## The third derivs
			pd2['sss']= Ixxx[i,j]
			pd2['dss']= Ixxy[i,j]
			pd2['dds']= Iyyx[i,j]
			pd2['ddd']= Iyyy[i,j]

			## Calculate the transfer functions
			y1 = data.ygroup[0].toTwoport(1,9)
			h1 = data.ygroup[0].transferFunction(1,9)
	   
			y2 = data.ygroup[1].toTwoport(1,9)
			h2 = data.ygroup[1].transferFunction(1,9)

			y3 = data.ygroup[2].toTwoport(1,9)
			h3 = data.ygroup[2].transferFunction(1,9)


			_v1,_v2,_v3 = [],[],[]
			_p0,_p1,_p2,_p3 = [],[],[],[]
			_i1,_i2,_i3 = [],[],[]
			
			for v in _v0:
			
				## Fundamental
				_i1.append(geti1(h1,pd2,v))
				_v1.append(h1*v)
				_p1.append(0.5*Rl*np.abs(_i1[-1])**2)
			 
				## Second Order Products
				_i2.append(geti2(h2,pd2,_v1[-1]))
				_v2.append(np.dot(la.inv(y2),[0,_i2[-1]])[0])
				_p2.append(0.5*Rl*np.abs(_i2[-1])**2)

				## Third Order Products
				ia=geti3(h3,pd2,_v1[-1])
				ib=2*getg2(h3,pd2)*_v1[-1]*_v2[-1]
				_i3.append(ia+ib)
				_v3.append(np.dot(la.inv(y3),[0,_i3[-1]])[0])
				_p3.append(0.5*Rl*np.abs(_i3[-1])**2)
				_p0.append((abs(v)*abs(v))/(8*Rs))
				
				
			dbm0 = 30+(10*np.log10(_p0))
			dbm1 = 30+(10*np.log10(_p1))
			dbm2 = 30+(10*np.log10(_p2))
			dbm3 = 30+(10*np.log10(_p3))

			ip2 = np.power(np.subtract(dbm1,dbm2),2)
			ip2 = list(ip2).index(min(ip2))
			_IP2[i,j] = dbm2[ip2]

			ip3 = np.power(np.subtract(dbm1,dbm3),2)
			ip3 = list(ip3).index(min(ip3))
			_IP3[i,j] = dbm3[ip3]
			_GH[i,j] = 10*np.log10(np.abs(h1*h1)*(Rs/Rl))
			_GP[i,j] = np.angle(h1, deg=True)
			
		## Plot IP2 and IP3 contours
		plotcontour(_vgs,_vds,_IP2,"$V_{gs}$ $(V)$", "$V_{ds}$ $(V)$",'$P_{IP2}$ $(dBm)$')
		plotcontour(_vgs,_vds,_IP3,"$V_{gs}$ $(V)$", "$V_{ds}$ $(V)$",'$P_{IP3}$ $(dBm)$')
		plotcontour(_vgs,_vds,_GH,"$V_{gs}$ $(V)$", "$V_{ds}$ $(V)$",'Power Gain $(dB)$')
		plotcontour(_vgs,_vds,_GP,"$V_{gs}$ $(V)$", "$V_{ds}$ $(V)$",'Voltage Phase Difference $(\phi)$')
		
		## Reset the circuit file
		Rs,Rl,Gm,Rds = 50.0,50.0,0.080348,228.653  
		data.setElement('Rs',Rs)
		data.setElement('Rl',Rl)
		data.setElement('Gm',Gm)
		data.setElement('Rds',Rds)
		plt.show()


	################################
	#	   IP2/IP3 vs Rs/Rl	   #
	################################
	if 1==0:
		vgs = np.arange(-3, 0.5, 0.01)
		vgd = np.arange(-3, 0.5, 0.01)
		_vgd,_vgs = np.meshgrid(vgd, vgs)
   
		## Get index of bias point (dumb way)
		VGD,VGS = -2.55,-0.55  
		for row,val in enumerate(vgd):
			if round(val,2)==VGD:
				break
		for col,val in enumerate(vgs):
			if round(val,2)==VGS:
				break

		######################################
		## CALCULATE CURRENT AND ITS DERIVS ##
		######################################

		## Get the current
		Ids = current(_vgd,_vgs)

		## Get first derivatives 
		Ix = _partialx(Ids,_vgd)
		Iy = _partialy(Ids,_vgs)

		## Get the second derivatives
		Ixx = _partialx(Ix,_vgd)
		Ixy = _partialy(Ix,_vgs)
		Iyy = _partialy(Iy,_vgs)

		## Get the third Derivatives
		Ixxx = _partialx(Ixx,_vgd)
		Ixxy = _partialy(Ixx,_vgs)
		Iyyx = _partialy(Ixy,_vgs)
		Iyyy = _partialy(Iyy,_vgs)

		## Create a dictionary of the partials
		## at the bias point.
		pd3 = {}
   
		## The constant term
		pd3[''] =Ids[col,row]

		## The first derivs
		pd3['d']  = Ix[col,row]
		pd3['s']  = Iy[col,row]

		## The second derivs
		pd3['dd'] = Ixx[col,row]
		pd3['ds'] = Ixy[col,row]
		pd3['ss'] = Iyy[col,row]

		## The third derivs
		pd3['ddd']= Ixxx[col,row]
		pd3['dds']= Ixxy[col,row]
		pd3['dss']= Iyyx[col,row]
		pd3['sss']= Iyyy[col,row]
		
		## Get a frequency list of the intermodulation products
		f1,f2,f3 = 1e10,1e10+1e5,1e10+2e5
		freqList = [f1,f2,f3]
		path=os.getcwd()+"/amplifier.cir"
		data = ms.frequencyAnalysis.fromFile(path, freqList)

		## Set Rs, Rl, Gm, Rds to those of the original amplifier
		## in case they have been modified by another simulation.
		Rs,Rl,Gm,Rds = 50.0,50.0,0.080348,228.653  
		data.setElement('Rs',Rs)
		data.setElement('Rl',Rl)
		data.setElement('Gm',Gm)
		data.setElement('Rds',Rds)

		npoints = 50
		_Rs = np.logspace(0,3,npoints)
		_Rl = np.logspace(0,3,npoints)

		_IP2,_IP3 =np.zeros([npoints,npoints]),np.zeros([npoints,npoints])
		_H,_P = np.zeros([npoints,npoints]),np.zeros([npoints,npoints])
		for i,Rs in enumerate(_Rs): 
			for j,Rl in enumerate(_Rl):

				## Set the source and load impedances
				data.setElement('Rs',Rs)
				data.setElement('Rl',Rl)
				
				print Rs, Rl

				## Perform a frequency analysis
				data = ms.frequencyAnalysis.fromFile(path, freqList)

				## Calculate the transfer functions
				y1 = data.ygroup[0].toTwoport(1,9)
				h1 = data.ygroup[0].transferFunction(1,9)
	   
				y2 = data.ygroup[1].toTwoport(1,9)
				h2 = data.ygroup[1].transferFunction(1,9)

				y3 = data.ygroup[2].toTwoport(1,9)
				h3 = data.ygroup[2].transferFunction(1,9)

				## calculate IP3
				_v0 = np.logspace(-5,2,300)

				_v1,_v2,_v3 = [],[],[]
				_p0,_p1,_p2,_p3 = [],[],[],[]
				_i1,_i2,_i3 = [],[],[]

				for v in _v0:
			
					## Fundamental
					_i1.append(geti1(h1,pd3,v))
					_v1.append(h1*v)
					_p1.append(0.5*Rl*np.abs(_i1[-1])**2)
					
					## Second Order Products
					_i2.append(geti2(h2,pd3,_v1[-1]))
					_v2.append(np.dot(la.inv(y2),[0,_i2[-1]])[0])
					_p2.append(0.5*Rl*np.abs(_i2[-1])**2)

					## Third Order Products
					ia=geti3(h3,pd3,_v1[-1])
					ib=2*getg2(h3,pd3)*_v1[-1]*_v2[-1]
					_i3.append(ia+ib)
					_v3.append(np.dot(la.inv(y3),[0,_i3[-1]])[0])
					_p3.append(0.5*Rl*np.abs(_i3[-1])**2)
					_p0.append(v*v/(8*Rs))

				dbm0 = 30+10*np.log10(_p0)
				dbm1 = 30+10*np.log10(_p1)
				dbm2 = 30+10*np.log10(_p2)
				dbm3 = 30+10*np.log10(_p3)

				ip2 = np.power(np.subtract(dbm1,dbm2),2)
				ip2 = list(ip2).index(min(ip2))
				_IP2[j,i] = dbm2[ip2]

				ip3 = np.power(np.subtract(dbm1,dbm3),2)
				ip3 = list(ip3).index(min(ip3))
				_IP3[j,i] = dbm3[ip3]
				_H[j,i] = 10*np.log10(np.abs(h1*h1)*(Rs/Rl))
				_P[j,i] = np.angle(h1, deg=True)

		plotcontour(_Rs,_Rl,_IP2,"$R_{s}$ $(\Omega)$", "$R_{l}$ $(\Omega)$",'$P_{IP2}$ $(dBm)$','log')
		plotcontour(_Rs,_Rl,_IP3,"$R_{s}$ $(\Omega)$", "$R_{l}$ $(\Omega)$",'$P_{IP3}$ $(dBm)$','log')
		plotcontour(_Rs,_Rl,_H,"$R_{s}$ $(\Omega)$", "$R_{l}$ $(\Omega)$",'Power Gain $(dB)$','lin')
		plotcontour(_Rs,_Rl,_P,"$R_{s}$ $(\Omega)$", "$R_{l}$ $(\Omega)$",'Voltage Phase Difference $(\phi)$','lin')


		## Reset the circuit file
		Rs,Rl,Gm,Rds = 50.0,50.0,0.080348,228.653  
		data.setElement('Rs',Rs)
		data.setElement('Rl',Rl)
		data.setElement('Gm',Gm)
		data.setElement('Rds',Rds)
		plt.show()

	################################
	#	 IP2/IP3 Static Values	#
	################################

	if 1==0:
		vgs = np.arange(-3, 0.5, 0.01)
		vgd = np.arange(-3, 0.5, 0.01)
		_vgd,_vgs = np.meshgrid(vgd, vgs)
   
		## Get index of bias point (dumb way)
		VGD,VGS = -2.50,-0.54  
		for row,val in enumerate(vgd):
			if round(val,2)==VGD:
				break
		for col,val in enumerate(vgs):
			if round(val,2)==VGS:
				break

		######################################
		## CALCULATE CURRENT AND ITS DERIVS ##
		######################################

		## Get the current
		Ids = current(_vgd,_vgs)

		## Get first derivatives 
		Ix = _partialx(Ids,_vgd)
		Iy = _partialy(Ids,_vgs)

		## Get the second derivatives
		Ixx = _partialx(Ix,_vgd)
		Ixy = _partialy(Ix,_vgs)
		Iyy = _partialy(Iy,_vgs)

		## Get the third Derivatives
		Ixxx = _partialx(Ixx,_vgd)
		Ixxy = _partialy(Ixx,_vgs)
		Iyyx = _partialy(Ixy,_vgs)
		Iyyy = _partialy(Iyy,_vgs)

		## Create a dictionary of the partials
		## at the bias point.
		pd3 = {}
   
		## The constant term
		pd3[''] =Ids[col,row]

		## The first derivs
		pd3['d']  = Ix[col,row]
		pd3['s']  = Iy[col,row]

		## The second derivs
		pd3['dd'] = Ixx[col,row]
		pd3['ds'] = Ixy[col,row]
		pd3['ss'] = Iyy[col,row]

		## The third derivs
		pd3['ddd']= Ixxx[col,row]
		pd3['dds']= Ixxy[col,row]
		pd3['dss']= Iyyx[col,row]
		pd3['sss']= Iyyy[col,row]
		
		## Get a frequency list of the intermodulation products
		f1,f2,f3 = 1e10,1e10+1e5,1e10+2e5
		freqList = [f1,f2,f3]
		path=os.getcwd()+"/amplifier.cir"
		data = ms.frequencyAnalysis.fromFile(path, freqList)

		## Set Rs, Rl, Gm, Rds to those of the original amplifier
		## in case they have been modified by another simulation.
		Rs,Rl,Gm,Rds = 50.0,50.0,0.080348,228.653  
		data.setElement('Rs',Rs)
		data.setElement('Rl',Rl)
		data.setElement('Gm',Gm)
		data.setElement('Rds',Rds)

		## Set the source and load impedances
		data.setElement('Rs',Rs)
		data.setElement('Rl',Rl)
	  
		## Perform a frequency analysis
		data = ms.frequencyAnalysis.fromFile(path, freqList)

		## Calculate the transfer functions
		y1 = data.ygroup[0].toTwoport(1,9)
		h1 = data.ygroup[0].transferFunction(1,9)
	   
		y2 = data.ygroup[1].toTwoport(1,9)
		h2 = data.ygroup[1].transferFunction(1,9)

		y3 = data.ygroup[2].toTwoport(1,9)
		h3 = data.ygroup[2].transferFunction(1,9)

		## calculate IP3
		_v0 = np.logspace(-5,2,300)

		_v1,_v2,_v3 = [],[],[]
		_p0,_p1,_p2,_p3 = [],[],[],[]
		_i1,_i2,_i3 = [],[],[]

		for v in _v0:
			
			## Fundamental
			_i1.append(geti1(h1,pd3,v))
			_v1.append(h1*v)
			_p1.append(0.5*Rl*np.abs(_i1[-1])**2)
					
			## Second Order Products
			_i2.append(geti2(h2,pd3,_v1[-1]))
			_v2.append(np.dot(la.inv(y2),[0,_i2[-1]])[0])
			_p2.append(0.5*Rl*np.abs(_i2[-1])**2)

			## Third Order Products
			ia=geti3(h3,pd3,_v1[-1])
			ib=2*getg2(h3,pd3)*_v1[-1]*_v2[-1]
			_i3.append(ia+ib)
			_v3.append(np.dot(la.inv(y3),[0,_i3[-1]])[0])
			_p3.append(0.5*Rl*np.abs(_i3[-1])**2)
			_p0.append(v*v/(8*Rs))

		dbm0 = 30+10*np.log10(_p0)
		dbm1 = 30+10*np.log10(_p1)
		dbm2 = 30+10*np.log10(_p2)
		dbm3 = 30+10*np.log10(_p3)

		ip2 = np.power(np.subtract(dbm1,dbm2),2)
		ip2 = list(ip2).index(min(ip2))
		_IP2 = dbm2[ip2]
		
		ip3 = np.power(np.subtract(dbm1,dbm3),2)
		ip3 = list(ip3).index(min(ip3))
		_IP3 = dbm3[ip3]

		print "IP2 = %f"%(_IP2)
		print "IP3 = %f"%(_IP3)
			   
		plt.figure(5)
		h1=plt.plot(dbm0,dbm1)
		h2=plt.plot(dbm0,dbm2)
		h3=plt.plot(dbm0,dbm3)
		plt.legend([h1,h2,h3],["$P_1$","$P_2$","$P_3$"],2)
		plt.xlabel("Input Power $(dBm)$")
		plt.ylabel("Output Power $(dBm)$")
		plt.show()

		## Reset the circuit file
		Rs,Rl,Gm,Rds = 50.0,50.0,0.080348,228.653  
		data.setElement('Rs',Rs)
		data.setElement('Rl',Rl)
		data.setElement('Gm',Gm)
		data.setElement('Rds',Rds)
		plt.show()
		