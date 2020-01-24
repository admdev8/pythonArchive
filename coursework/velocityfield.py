#!/usr/bin/env python
import numpy as np
import MonteCarlo as MC
import cPickle as pickle
import multiprocessing as mp
import matplotlib.pyplot as plt

## For Plotting
from matplotlib import rc
rc("font",**{"family":"sans−serif","sans−serif":["Helvetica"],"size":20})
rc("text", usetex=True)
rc("xtick", labelsize=20)
rc("ytick", labelsize=20)
rc("legend", fontsize=20)
params = {"legend.fontsize": 16}
plt.rcParams.update(params) 

################################
# MULTIPROCESS CODE 
#
class AsyncFactory:

	def __init__(self, func, cb_func):

		self.func, self.cb_func, self.pool = func, cb_func, mp.Pool()

	def call(self,*args, **kwargs):

		self.pool.apply_async(self.func, args, kwargs, self.cb_func)

	def wait(self):

		self.pool.close()
		self.pool.join()

	def run_MonteCarlo(F, NPOINTS, mode, EMAX = 1.0):

		E = np.linspace(0.0,EMAX,100)
		M = MC.MonteCarlo(E, F, NPOINTS, mode)
		ELECTRON = M.run()
		ELECTRON.buildData()
		print F ,−1*ELECTRON.VMEAN
		return (ELECTRON.VMEAN ,F, ELECTRON.VALLEY)
		result_list = []

	def log_result(result):

		result_list.append(result)


###############################
# Running the Simulator 
#
if __name__ == "__main__":

	if False:

		FIELD,NPOINTS = [1000, 3000, 10000, 30000], 100000
		async_MonteCarlo = AsyncFactory(run_MonteCarlo,log_result)

		for F in FIELD: 

			async_MonteCarlo.call(F, NPOINTS, "FULL", 1.0)
			async_MonteCarlo.wait()
			result_list = sorted(result_list,key=lambda x: x[1])
			VMEAN = [abs(E[0]) for E in result_list]
			VELOCITY = [E[2] for E in result_list]
		
		# Dump the results so they can be plotted later
		data = {"FIELD":FIELD,"VELOCITY":VELOCITY}
		output = open("100k_full2.pkl", "wb")
		pickle.dump(data, output)
		plt.plot(VELOCITY[0])
		plt.show()
		

	if False:
		
		FIELD,NPOINTS = np.linspace(10000.0, 30000.0, 50), 100000
		async_MonteCarlo = AsyncFactory(run_MonteCarlo,log_result)

		for F in FIELD: 

			async_MonteCarlo.call(F, NPOINTS, "AVG", 1.0)
			 async_MonteCarlo.wait()
			result_list = sorted(result_list,key=lambda x: x[1])

			# result_list is a list of electron objects. We would
			# like to get out the mean velocity of each electron
			# and plot it as a function of electric field
			VMEAN = [abs(E[0]) for E in result_list]

			# Also get the valley occupancies
			VALLEY_G = [E[2]["G"] for E in result_list]
			VALLEY_L = [E[2]["L"] for E in result_list]

		# Dump the results so they can be plotted later
		data = {"FIELD":FIELD,"VMEAN":VMEAN,"VALLEY_G":VALLEY_G,"VALLEY_L":VALLEY_L}
		output = open("30kV_100k_avg.pkl", "wb")
		pickle.dump(data, output)

		# Plot the data when done
		plt.figure(1)
		plt.plot(FIELD/1000, VMEAN)
		plt.xlabel("Electric Field $(kV/cm)$")
		plt.ylabel("Electron Velocity $(cm/s)$")

		plt.figure(2)
		h1=plt.plot(FIELD/1000, VALLEY_G)
		h2=plt.plot(FIELD/1000, VALLEY_L)
		plt.xlabel("Electric Field $(kV/cm)$")
		plt.ylabel("Valley Occupancy")
		plt.legend([h1,h2], ["$\Gamma$", "L"])
		plt.show()
