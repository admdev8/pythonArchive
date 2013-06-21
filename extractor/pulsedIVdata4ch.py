#!/usr/bin/env python 

#import libUtils
#import physics

# For data extraction 
import matplotlib.pyplot as plt
import numpy as np

# For pickling
import cPickle as pickle
from collections import defaultdict

# For getting files from root directory and printing errors
from os import chdir
from os import listdir
from os.path import isfile, join
import operator
import sys

from matplotlib import rc
rc('xtick', labelsize=20) 
rc('ytick', labelsize=20) 
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':20})
plt.rcParams.update({'legend.fontsize':16,'legend.linewidth':2})
rc('text', usetex=True)

# Helper method
def zigzag(seq): 
  return seq[::2], seq[1::2]


def filter_value(_list, value):
    for x,y,z in _list:
        if x == value:
            yield y,z

# Pulse Ladder Class: This contains the core extraction 
# Algorithm for a two channel pulse ladder
class pulseLadder(object): 

    def __init__(self, _file, _pulseLadder, _voltage, _current): 
        #self.name = _file
        self.pulseLadder = _pulseLadder
        #self.voltage = _voltage
        #self.current = _current

        # non-linear least squares fit
        #fit_guess = [.01,0.05,.02]
        #fit_step = [.0001,.0001,.0001]
        #convergence,maxstep = 5e-6, 10000

        #self.vfit,self.ifit,self.resistance,self.p_found = physics.ivTanhFit(self.voltage,self.current,fit_guess,fit_step,convergence,maxstep)

        # Checking for resistance curves
        #plt.figure(1)
        #plt.plot(self.voltage,self.current)
        #plt.plot(self.vfit,self.ifit)
        #plt.figure(2)
        #plt.plot(self.vfit,self.resistance)
        #plt.show()

    @classmethod
    def datFromFile(cls, _file): 
        name = _file
        if ".pkl" in _file: 
            IVinput = open(_file,'rb')
            _pulseLadder = pickle.load(IVinput)
            IVinput.close()

        else: 
            sys.stderr.write("No .pkl file, check directory") 
            exit(1)

        # Obtain Voltages: Each key is a Voltage (V*V)
        #_voltage = []
        #for i,key in enumerate(_pulseLadder): 
        #    _vConvolve = list(np.convolve(_pulseLadder[key][1],_pulseLadder[key][1],mode='full'))
        #    _vConvolve,e = zigzag(_vConvolve)
        #    _voltage.append(_pulseLadder[key][1][_vConvolve.index(max(_vConvolve))])
        # Obtain Currents: Each key is a Voltage (V*I)

        _current,_vgreal, _vdreal, _vg, _vd = [],[],[]
        for i,key in enumerate(_pulseLadder): 

            # Perform convolution to find pulses
            _iConvolve = list(np.convolve(_pulseLadder[key][2],_pulseLadder[key][1],mode='full'))
            _iConvolve,e = zigzag(_iConvolve)
            _current.append(_pulseLadder[key][2][_iConvolve.index(max(_iConvolve))])
            # Use same maximum to get the real  
            _vdreal.append(_pulseLadder[key][1][_iConvolve.index(max(_iConvolve))])
            _vgreal.append(_pulseLadder[key][3][_iConvolve.index(max(_iConvolve))])
            # Append the setpoint values
            _vg.append(float(key.split()[0]))
            _vd.append(float(key.split()[1]))
                 
        data = sorted(zip(_vg,_vd,_current),key=operator.itemgetter(0))
        vz,hz = [],[]
        for v in sorted(set(_vg)):
             result= list(filter_value( data, v))
             x = [list(i) for i in zip(*sorted(result,key=operator.itemgetter(0)))]
             h1 = plt.plot(x[0],x[1], 'o',color = plt.get_cmap('jet')(float(v)/(len(sorted(set(_vg)))-1)))
             vz.append("$V_g = %sV$"%v)
             hz.append(h1)

        plt.legend(hz,vz,loc=2)
        plt.show()

        #return cls(_file, _pulseLadder, _voltage, _current)

    # For legend include plt.legend(hlist,self.pulseLadder.keys())   
    def plotPulseLadders(self, title): 

        # Plot Pulses
        plt.figure(1)
        hlist = []
        for i,key in enumerate(self.pulseLadder): 
            h1 = plt.plot(np.multiply(1000000,self.pulseLadder[key][0]),self.pulseLadder[key][1])
            hlist.append(h1)
            plt.title(title)
            plt.xlabel(r'Time $(\mu s)$')
            plt.ylabel(r'Voltage $(V)$')
            
        plt.figure(2)
        hlist = []
        for i,key in enumerate(self.pulseLadder): 
            h1 = plt.plot(np.multiply(1000000,self.pulseLadder[key][0]),np.multiply(1000,self.pulseLadder[key][2]))
            hlist.append(h1)
            plt.xlabel(r'Time $(\mu s)$')
            plt.ylabel(r'Current $(mA)$')
            plt.xlim(3,8)
        plt.show()

    def plotConvolveLadders(self, title):
        
        # Perform Convolution to find peaks 
        plt.figure(3)
        hlist = []
        for i,key in enumerate(self.pulseLadder): 
            Vfind = list(np.convolve(self.pulseLadder[key][1],self.pulseLadder[key][1],mode='full'))
            Vfind,e = zigzag(Vfind)
            h1 = plt.plot(self.pulseLadder[key][0],Vfind)
            hlist.append(h1)
            plt.title(str(title))
            plt.xlabel("Time (s)")
            plt.ylabel("V*V")
        
        plt.figure(4)
        hlist = []
        for i,key in enumerate(self.pulseLadder): 
            Ifind = list(np.convolve(self.pulseLadder[key][2],self.pulseLadder[key][1],mode='full'))
            Ifind,e = zigzag(Ifind)
            h1 = plt.plot(self.pulseLadder[key][0],Ifind)
            hlist.append(h1)
            plt.title(str(title))
            plt.xlabel("Time (s)")
            plt.ylabel("V*I")
        plt.show()

    def plotIV(self, title):
        h1 = plt.plot(self.voltage, np.multiply(1000,self.current),'bo')
        #h2 = plt.plot(self.vfit, np.multiply(1000,self.ifit),'k')
        plt.xlabel(r'Voltage $(V)$')
        plt.ylabel(r'Current $(mA)$')
        #plt.legend([h1,h2], [r'IV data',r'$I_0 [1 + \lambda V] \tanh(\beta V)$'],loc=2, scatterpoints = 1, numpoints = 1)
        plt.show()

        

# This is simply a list of pulse ladder objects
# All it can do is store the IV curves in a dict 
# and plot them if needed. It does not do much
class pulseList(object): 
    def __init__(self,_pulseList): 
        self.pulseList = _pulseList
    
    def __getitem__(self,key): 
        return self.pulseList[key]

    @classmethod
    def listFromDir(cls, _rootpath):
        ext = '.dat'
        _pulseList = {}
               
        # Get file list in rootpath
        files = [f for f in listdir(_rootpath) if isfile(join(_rootpath,f)) ]
 
        # enter rootpath, open each file, and save iv
        chdir(_rootpath)
        for _file in files:
            print "processing "+_file
            _pulseLadder = pulseLadder.datFromFile(_file)
            _pulseList[_file] = _pulseLadder

        return cls(_pulseList) 

    def plotPulseListIV(self, title): 
        hlist = []
        for i,key in enumerate(self.pulseList): 
            plt.plot(self.pulseList[key].vfit,np.multiply(1000,self.pulseList[key].ifit),'k')
            h1 = plt.plot(self.pulseList[key].voltage,np.multiply(1000,self.pulseList[key].current),'o')
            hlist.append(h1)
        
        plt.title(r'H-intercalated IV 300K', fontsize=16)
        plt.xlabel(r'Voltage $(V)$',fontsize=16)
        plt.ylabel(r'Current $(mA)$',fontsize=16)
        plt.legend(hlist,self.pulseList.keys(),loc=2, scatterpoints = 1, numpoints = 1)   
        plt.ylim([0,50])
        plt.show()


#########################################################################################################
if __name__ == "__main__": 
    #########################
    # Single File Operation #
    #########################
    _file = "../data/test/testdata4ch_4.pkl"
    _300K = pulseLadder.datFromFile(_file)
    






    #_300K.plotPulseLadders(r'')   
    #_300K.plotConvolveLadders("Pulse 300K")
    #_300K.plotIV(r'')

    ############################
    # From Directory Operation #
    ############################
    #_rootpath = "./wunjo/K/dat"
    #_D = pulseList.listFromDir(_rootpath)
    #_D.plotPulseListIV("as-grown")
    
    #_data = defaultdict(list)
    #for key in _D.pulseList.keys():
    #  _data[key].append(_D.pulseList[key].vfit)
    #  _data[key].append(_D.pulseList[key].ifit)
    #  _data[key].append(_D.pulseList[key].resistance)
     

    #_pkl = open("./pkl/K_wunjo.pkl",'wb')
    #pickle.dump(_data,_pkl)
    


   
