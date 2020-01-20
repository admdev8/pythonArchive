#!/usr/bin/env python 
import os, sys, inspect

os.chdir('../drivers')
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

import pulsedIV_driver as pd
import matplotlib.pyplot as plt 
import numpy as np
import cPickle as pickle
import time
import math

class pulsedIV_2ch_meas(object): 

    def __init__(self, sourceGPIB, scopeGPIB, pulseGPIB): 

        self.source = pd.source_66xx(sourceGPIB)
        self.scope = pd.scope_tds7104(scopeGPIB)
        self.pulse = pd.tabor_8500(pulseGPIB)

    def setPulses(self, length, delay): 
        self.length = length
        self.delay  = delay
        self.pulse.setPulseA(self.length, self.delay)
         
    def measure(self, vList, numavg): 

        self.data = {}
        # Pulse through the voltage list
        for v in vList: 

            print "Bias Settle:",round(v,1)
            self.source.setVoltage(0,round(v,1),0.5) 
            self.scope.setAverage(numavg)
            
            # Stabilize Voltage
            time.sleep(1)
           
            # Send pulses
            self.pulse.setVoltageA(0.0,2.0)
            self.pulse.triggerA()
            
            # Buffer sleep
            time.sleep(math.ceil(numavg/100)+1)
            self.pulse.offPulseA()
            self.scope.setAcqOff()

            # Scope Data
            tch1, vch1 = self.scope.chx_data(1)
            tch2, vch2 = self.scope.chx_data(2)

            self.data["%s"%v] = [tch1,vch1,vch2]

class pulsedIV_4ch_meas(object): 

    def __init__(self, sourceGPIB, scopeGPIB, pulseGPIB): 

        self.source = pd.source_66xx(sourceGPIB)
        self.scope = pd.scope_tds7104(scopeGPIB)
        self.pulse = pd.tabor_8500(pulseGPIB)

    def setPulses(self, length, delay): 
        self.length = length
        self.delay  = delay
        self.pulse.setPulseA(self.length, self.delay)
        
    def measure(self, vList, gList, numavg): 

        self.data = {}
        # Pulse through the voltage lists
        for g in gList: 
            print "Vg = %s"%round(g,1)
            print "--------------------"

            self.pulse.setVoltageA(0,round(g,1))
            self.pulse.triggerA()
            
            for v in vList: 

                print "Vd = %s"%round(v,1)
                self.source.setVoltage(0,round(v,1),0.5) 
                self.scope.setAverage(numavg)
                # Stabilize Voltage
                time.sleep(1)
            
                # Buffer sleep
                time.sleep(math.ceil(numavg/100)+1)
                
                self.scope.setAcqOff()

                # Scope Data
                tch1, vch1 = self.scope.chx_data(1)
                tch2, vch2 = self.scope.chx_data(2)
                tch1, vch3 = self.scope.chx_data(3)
                tch2, vch4 = self.scope.chx_data(4)
                self.data["%s %s"%(g,v)] = [tch1,vch1,vch2,vch3,vch4]
            
            self.pulse.offPulseA()
            print "\n"
    
if __name__ == "__main__": 

    if 1 == 1: 
        root = "/home/hadit/Desktop/chipA/PulsedIV/test/"
        path = "test.pkl"

        print "Saving data in:%s"%root+path 

        vList = np.linspace(0,9,31)
        # Specify GPIB addresses. 
        meas = pulsedIV_2ch_meas(6,15,9)
        # Specify pulse length (us) and delay (ms)
        meas.setPulses(1,1)
        meas.scope.setScale(1,10)
        meas.scope.setScale(2,0.01)
        # Set number of averages and perform measurement
        meas.measure(vList,2048)


        plt.figure(1)
        for key, value in meas.data.iteritems():
            plt.plot(value[0],value[1])
        plt.title("Test Pulses after USB Atrocity")
        plt.xlabel("Time (s)")
        plt.ylabel("Voltage (V)")
   
        plt.figure(2)
        for key, value in meas.data.iteritems():
            plt.plot(value[0],value[2])
        plt.title("Test Pulses after USB Atrocity")
        plt.xlabel("Time (s)")
        plt.ylabel("Current (A)")
        plt.show()

        pickle.dump(meas.data,open(root+path, "wb"))

     ## FOUR CHANNEL MEASUREMENTS ##
    if 1 == 0: 
        vList,gList = np.linspace(0,5,26),np.linspace(0,4,9)
        meas = pulsedIV_4ch_meas(6,15,9) 
        meas.setPulses(2,1)
        meas.scope.setScale(2,0.05)
        meas.measure(vList,gList,512)

        plt.figure(1)
        for key, value in meas.data.iteritems():
            plt.plot(value[0],value[1])
        plt.title("Drain Voltage")
        plt.xlabel("Time (s)")
        plt.ylabel("Voltage (V)")
   
        plt.figure(2)
        for key, value in meas.data.iteritems():
            plt.plot(value[0],value[2])
        plt.title("Drain Current")
        plt.xlabel("Time (s)")
        plt.ylabel("Current (A)")
        plt.show()

        plt.figure(3)
        for key, value in meas.data.iteritems():
            plt.plot(value[0],value[3])
        plt.title("Gate Voltage")
        plt.xlabel("Time (s)")
        plt.ylabel("Voltage (V)")
   
        plt.figure(4)
        for key, value in meas.data.iteritems():
            plt.plot(value[0],value[4])
        plt.title("Gate Current")
        plt.xlabel("Time (s)")
        plt.ylabel("Current (A)")
        plt.show()

        pickle.dump(meas.data,open("../data/test/testdata4ch_4.pkl", "wb"))


        
    



