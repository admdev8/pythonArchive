#!/usr/bin/env python 

# Plotting Libs
import matplotlib.pyplot as plt  
import cPickle as pickle
import pprint
import numpy as np
import os

# Plotting font management
from matplotlib import rc
rc('xtick', labelsize=20) 
rc('ytick', labelsize=20) 
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':18})
plt.rcParams.update({'legend.fontsize':16,'legend.linewidth':2})
rc('text', usetex=True)

class ramanMap(object): 

    def __init__(self, pathname, size):

        self.size = size

        pkl_file = open(pathname, 'rb')
        mapData = pickle.load(pkl_file)

        keys,cent,line,ampl,nois = [],[],[],[],[]
        for key in mapData.keys(): 
            keys.append(int(key)) 
            cent.append(mapData[key][0])
            line.append(mapData[key][1])
            ampl.append(mapData[key][2])
            nois.append(mapData[key][3])
            
        self.peakcenter = self.zipMap(keys,cent)
        self.linewidth = self.zipMap(keys,line)
        self.amplitude = self.zipMap(keys,ampl)
        self.fitnoise = self.zipMap(keys,nois)
        pkl_file.close()

    def zipMap(self, keys, vals): 
        return [list(x) for x in zip(*sorted(zip(keys, vals), key=lambda pair: pair[0]))]

    def plotMap(self, data): 

        step = 0.3
        x,y = np.zeros((self.size,self.size)),np.zeros((self.size,self.size))
    
        for i in range(self.size):
            for j in range(self.size): 
                index = i*self.size+j
                try: 
                    x[i][j] = index
                    y[i][j] = data[1][index]
                except: 
                     y[i][j] = np.nan
    
        plt.imshow(y, origin='lower', interpolation='nearest',cmap = plt.get_cmap('RdBu'), extent=[0,self.size*step,0,self.size*step])
        plt.xlabel(r'x-position $($\mu m$)$')
        plt.ylabel(r'y-position $($\mu m$)$')
        plt.colorbar()
        
    

pathname = "./maps/BML39I/2Dpeak_third.pkl" 
mapLIU   = ramanMap(pathname, 15)
plt.figure(1)
plt.title(r'2D Peak Position Map $(cm^{-1})$')
mapLIU.plotMap(mapLIU.peakcenter)
plt.figure(2)
plt.title(r'2D Peak Linewidth Map $(cm^{-1})$')
mapLIU.plotMap(mapLIU.linewidth)
plt.figure(3)
plt.title(r'2D Peak Lorentz Fit Noise')
mapLIU.plotMap(mapLIU.fitnoise)
plt.show()

pathname = "./maps/BML39I/Gpeak_LIU.pkl" 
mapLIUG   = ramanMap(pathname, 15)
#plt.figure(4)
#mapLIUG.plotMap(mapLIUG.linewidth)
#plt.figure(5)
#mapLIUG.plotMap(mapLIUG.peakcenter)
#plt.figure(6)
#mapLIUG.plotMap(mapLIUG.fitnoise)


pathname = "./maps/poland/poland.pkl" 
poland = ramanMap(pathname, 15)
pathname = "./maps/poland/Gpeak_poland.pkl" 
polandG = ramanMap(pathname, 15)


plt.figure(7)
#plt.title(r'2D Peak Position vs. Linewidth')
h1 = plt.plot(mapLIU.peakcenter[1], mapLIU.linewidth[1],'o')
h2 = plt.plot(poland.peakcenter[1], poland.linewidth[1],'ro')
plt.xlabel(r'2D Peak Position $(cm^{-1})$')
plt.ylabel(r'2D Peak Linewidth $(cm^{-1})$')
plt.legend([h1,h2],['Sample 1', 'Sample 2'], loc = 2,scatterpoints = 1, numpoints = 1)
plt.show()
#plt.figure(8)
#plt.title(r'G peak Position vs. Linewidth')
#h1 = plt.plot(mapLIUG.peakcenter[1], mapLIUG.linewidth[1],'o')
#h2 = plt.plot(polandG.peakcenter[1], polandG.linewidth[1],'ro')
#plt.xlabel(r'G Peak Position $(cm^{-1})$')
#plt.ylabel(r'G Peak Linewidth $(cm^{-1})$')
#plt.legend([h1,h2],['LIU', 'ITME'], loc = 2,scatterpoints = 1, numpoints = 1)
#plt.show()
