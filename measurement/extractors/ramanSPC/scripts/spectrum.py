#!/usr/bin/env python 
import cPickle as pickle
import copy
import math
import os
import re

# Import ndfit
import ndfit as ndf

# Numerics and plotting
import matplotlib.pyplot as plt
import numpy as np

# This script contains methods to extract and manipulate Raman spectroscopy ascii data
# Note that spectral data is typucally contained in .prn files.

# Class to hold a spectrum
class spectrum: 
    
    def __init__(self, l, c, substrate = None): 
     
        self.laser = 5319.8
        self.wavelength, self.counts = l,c
        self.r_shift = [1e8*((1./self.laser)- (1./i)) for i in self.wavelength]
        self.substrate_raw = substrate[1]

    @classmethod
    def fromFile(cls, _file, substrate = None): 

        l,c = [],[]
        tmp = open(_file, "rb") 
        for line in tmp: 
            l.append(float(line.split()[0])) 
            c.append(float(line.split()[1]))
  
        if substrate is not None: 
            pkl_file = open(substrate, 'rb')
            substrate = pickle.load(pkl_file)
            spectrum = cls(l,c, substrate)
            return spectrum

        else:
            spectrum = cls(l,c, substrate)
            return spectrum

    def despike(self): 

        waveList = [5,10,15,20]
        self.counts_despiked = copy.deepcopy(self.counts)

        for val in waveList:
            _wx = np.linspace(-1,1,val)
            _wy = 2*np.sinc(2*_wx)-np.sinc(_wx)
    
            y = np.convolve(_wy/_wy.sum(),self.counts,mode='same')
            window = 3
            for i,val in enumerate(np.sign(self.counts*y)): 
                if val == -1: 
                    if i < len(self.counts_despiked)-window:
                        self.counts_despiked[i] = (self.counts_despiked[i+window]+self.counts_despiked[i-window])/2

    
        # A smoothed and despiked version
        self.counts_smoothed_despiked = self.smooth(self.counts_despiked)

    # Method to smooth data
    def smooth(self, data=None): 
        window_len, beta = 30, 30
        w = np.kaiser(window_len,beta)
        
        if data is None: 
            return np.convolve(w/w.sum(),self.counts,mode='same')
        else:
            return np.convolve(w/w.sum(),data,mode='same')
 
    # Method to select a window of the raman data for operations
    def selectWindow(self, x1, x2): 
        for i in self.r_shift: 
            if i > x1:
                index1 = self.r_shift.index(i) 
                break
                
        for i in self.r_shift: 
            if i > x2:
                index2 = self.r_shift.index(i)
                break

        self.r_window = self.r_shift[index1:index2]
        self.c_window = self.counts[index1:index2] 
        self.c_window_smoothed_despiked = self.counts_smoothed_despiked[index1:index2]
        
        if self.substrate_raw is not None:
            self.c_window_smoothed_despiked_nobg = self.counts_smoothed_despiked_nobg[index1:index2]
        
    def fitLorentzian(self,fit_guess,fit_step,convergence,maxdepth): 
                
        # Store Data
        if self.substrate_raw is not None:
            _data = [self.r_window, self.c_window_smoothed_despiked_nobg]
        else: 
            _data = [self.r_window, self.c_window_smoothed_despiked]
 
        # Hyperbolic model
        fitfunc = lambda p,d: ((p[2]/math.pi)*(p[1]))/((d[0]-p[0])**2 + p[1]**2) 
        errfunc = lambda p,d: fitfunc(p,d)-d[1]
  
        # Run recursive fit with lattice throttling
        fit = ndf.ndfit(_data, fitfunc, errfunc, fit_guess, fit_step, convergence, maxdepth)
        fit.lattice_select("fco")
        #fit.fit_throttling(1e-3, 2, 1e-3, 2)
        fit.fit_run()
        
        # Print recursion informtiton 
        print "recursion depth:"+str(fit.depth) 
        print "convergence:"+str(fit.res_vec[-1])
        print "noise amplitude:"+str(fit.res_vec[-1]/len(_data[0])) 
        print "params:"+str(fit.p_found) 
        
        return fit.p_found, fit.res_vec[-1]/len(_data[0]) 

        # Generate the nd_curve (the result)
        #curve = fit.nd_curve()
        #plt.figure(1)
        #plt.plot(curve[0],curve[1])
        #plt.plot(self.r_window, self.c_window_smoothed)
        #plt.show()

    def subtractBackground(self, xmin=None,xmax=None): 

        if xmin is None: xmin = 0
        if xmax is None: xmax = len(self.substrate_raw)

        scale, sum_vec = np.linspace(0.5,1.5,1000),[1e13]
        for s in scale:
            _tmp = [s*i for i in self.substrate_raw]
            _sum = sum([(i-j)**2 for i,j in  zip(self.counts_smoothed_despiked[xmin:xmax],_tmp[xmin:xmax])])
            sum_vec.append(_sum)
            if _sum > sum_vec[-2]: 
                break
            
        self.substrate_scaled = _tmp
        self.counts_smoothed_despiked_nobg = self.smooth([(i-j) for i,j in  zip(self.counts_smoothed_despiked,_tmp)])

    def plotSpectrum(self,xmin,xmax): 
        plt.figure(1) 
        plt.subplot(311)
        plt.title(r'Raman Spectra:H-Intercalated Graphene')
        plt.ylabel(r'Raw Data')
        plt.plot(self.r_shift, self.counts)
        plt.xlim([xmin,xmax])
        plt.subplot(312)
        plt.ylabel(r'Wavelet Despiked')
        plt.plot(self.r_shift, self.counts_despiked)
        plt.xlim([xmin,xmax])
        plt.subplot(313)
        plt.ylabel(r'Hanning Smoothed')
        plt.xlabel(r'Raman Shift $(cm^{-1})$')
        plt.plot(self.r_shift, self.counts_smoothed_despiked)
        plt.xlim([xmin,xmax])

    def plotSmoothed(self, xmin, xmax): 
        plt.figure(1) 
        plt.ylabel(r'Intensity $(counts)$')
        plt.xlabel(r'Raman Shift $(cm^{-1})$')
        h = plt.plot(self.r_shift, self.counts_smoothed_despiked, 'k')
        plt.xlim([xmin,xmax])
        return h

    def plotSubstrate(self, xmin, xmax): 
        plt.figure(1) 
        plt.ylabel(r'Intensity $(counts)$')
        plt.xlabel(r'Raman Shift $(cm^{-1})$')
        h = plt.plot(self.r_shift, self.substrate_raw,'0.65')
        plt.xlim([xmin,xmax])
        return h

    def plotSubstrateScaled(self, xmin, xmax): 
        plt.figure(1) 
        plt.ylabel(r'Intensity $(counts)$')
        plt.xlabel(r'Raman Shift $(cm^{-1})$')
        h = plt.plot(self.r_shift, self.substrate_scaled,'0.8')
        plt.xlim([xmin,xmax])
        return h
    
    def plotSmoothedNoSub(self, xmin, xmax): 
        plt.figure(1) 
        plt.ylabel(r'Intensity $(counts)$')
        plt.xlabel(r'Raman Shift $(cm^{-1})$')
        h = plt.plot(self.r_shift, self.counts_smoothed_despiked_nobg, 'r')
        plt.xlim([xmin,xmax])
        return h

# Spectral group object
class spectralGroup: 

    def __init__(self, d): 
        self.spectralGroup = d

    @classmethod
    def readDir(cls, pathname, substrate=None, val=0): 
        
        dataList = []
        # Get a list of files for the  map
        filelist = [f for f in os.listdir(pathname) if os.path.isfile(os.path.join(pathname,f))]
        
        # Get the mapidentifier for each file
        for f in filelist: 
            spec = [re.findall('\d+',f)[val],spectrum.fromFile(os.path.join(pathname,f),substrate)]
            dataList.append(spec)

        # Initialize the class
        spectralGroup = cls(dataList)
        return spectralGroup

    def despike(self):
        for f in self.spectralGroup: 
            f[1].despike()

    def subtractBackground(self,xmin,xmax):
        for f in self.spectralGroup: 
            f[1].subtractBackground(xmin,xmax)

    def averageSpectra(self): 
        
        avg = [0]*len(self.spectralGroup[0][1].r_shift)
        for f in self.spectralGroup: 
            avg = [i+j for i,j in zip(avg,f[1].counts_smoothed_despiked)]

        avg = [i/len(self.spectralGroup) for i in avg]
        return self.spectralGroup[0][1].r_shift,avg
        
    def selectWindow(self, x1, x2):
        for f in self.spectralGroup: 
            f[1].selectWindow(x1,x2)

    
# Static method to generate background data
def generateAverage(dirname, _id = 0):
      x  = spectralGroup.readDir(dirname,_id)
      x.despike()

      substrate = [[],[]]
      substrate[0],substrate[1] = x.averageSpectra()

      pathname = os.getcwd()
      pklpath = os.path.join(pathname,'substrate.pkl')
      output = open(pklpath, 'wb')
      #Pickle the dict using the highest protocol available.
      pickle.dump(substrate, output, -1)
      output.close()
      
      return pklpath

if __name__ == "__main__": 


    # Method to generate substrate average
    # substratePath = "./substrate/prn" 
    # pklpath = generateAverage(substratePath,1)
   
    if 1 == 0:

        substrate = "./substrate/substrate.pkl"
        filename = "./BML39I/prn/BML39I-rmap100.prn"
        x = spectrum.fromFile(filename, substrate)
        x.despike()
        plt.subplot(212)
        #plt.title(r'Grapehene Spectrum Isolation')
        x.subtractBackground(200,300)
        h1 = x.plotSmoothed(1570,1850)
        h2 = x.plotSubstrate(1570,1850)
        h4 = x.plotSmoothedNoSub(1570,1850)
        h3 = x.plotSubstrateScaled(1570,1850)
        plt.subplot(211)
        h1 = x.plotSmoothed(1570,2850)
        h2 = x.plotSubstrate(1570,2850)
        h4 = x.plotSmoothedNoSub(1570,2850)
        h3 = x.plotSubstrateScaled(1570,2850)
        plt.xlabel("")
        plt.legend([h1,h2,h3,h4],['Despiked Spectrum','Despiked SiC Background','Scaled SiC', 'Graphene'], loc = 9)
        plt.show()


    if 1 == 1:

        substrate = "./substrate/substrate.pkl"
        filename = "./BML39I/prn"

        x = spectralGroup.readDir(filename,substrate,1)
        x.despike()
        x.subtractBackground(200,300)
        x.selectWindow(1570,1620)
      
        _2Dpeak = {}
        #Stored as center, linewidth, norm
        fitGuess,fitStep,convergence,maxdepth = [1594,1,1100],[.1,.1,2],5e5,1024
        for f in x.spectralGroup: 
            params, noise = f[1].fitLorentzian(fitGuess,fitStep,convergence,maxdepth)
            _2Dpeak[f[0]] = [params[0],2*params[1], params[2], noise]
                
        output = open('./Gpeak_LIU.pkl', 'wb')
        #Pickle the dict using the highest protocol available.
        pickle.dump(_2Dpeak, output, -1)
        output.close()

    if 1 == 0:

        filename = "./BML39I/prn/BML39I-rmap100.prn"
        x = spectrum.fromFile(filename)
        x.despike()
        x.subtractBackground()
        x.plotSpectrum()

        filename = "./BML39I/prn/BML39I-rmap79.prn"
        y = spectrum.fromFile(filename)
        y.despike()
        y.plotSpectrum()

        filename = "./BML39I/prn/BML39I-rmap87.prn"
        z = spectrum.fromFile(filename)
        z.despike()
        z.plotSpectrum()
        plt.show()


    if 1 == 0:

    # Operation on a window with lorentzian fitting
        x.selectWindow(1570,1620)
        y.selectWindow(1570,1620)
        z.selectWindow(1570,1620)
       
        plt.figure(3)
        plt.title("Graphene G Peak Intensity")
        h1=plt.plot(x.r_window, x.c_window_smoothed_despiked)
        h2=plt.plot(y.r_window, y.c_window_smoothed_despiked)
        h3=plt.plot(z.r_window, z.c_window_smoothed_despiked)
        plt.xlabel(r'Raman Shift $(cm^{-1})$')
        plt.ylabel(r'Intensity $(counts)$')
        plt.legend([h1,h3,h2],["1 Layer", "1-2 Layers", "2 Layers"])
        plt.show()


    if 1 == 0:
        filename = "./poland/prn/"
        x = spectralGroup.readDir(filename,0)
        x.despike()
        x.selectWindow(2600,2800)
  
    
        _2Dpeak = {}
        #Stored as center, linewidth, norm
        fitGuess,fitStep,convergence,maxdepth = [2700,15,3000],[1,.1,10],5e5,1024
        for f in x.spectralGroup: 
            params, noise = f[1].fitLorentzian(fitGuess,fitStep,convergence,maxdepth)
            _2Dpeak[f[0]] = [params[0],2*params[1], params[2], noise]
        
        
        output = open('./poland.pkl', 'wb')
        #Pickle the dict using the highest protocol available.
        pickle.dump(_2Dpeak, output, -1)
        output.close()
