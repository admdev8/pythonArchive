#!/usr/bin/env python 
import math
import visa
import re
import os 
import time
import numpy as np
import matplotlib.pyplot as plt
import ndfit


class hp_4145B:

    def __init__(self, GPIB): 
        ## Initialize the dude
        self.analyzer = visa.instrument("GPIB::%s"%GPIB)


    def channelDefinition(self):
 
        ## Select channel [DE]finition page
        self.analyzer.write("DE")
       
        ## Perform Channel Definition Steps. We need 
        ## Two channels for a basic DC measurement.
        ## 1 = Voltage Mode
        ## 2 = Current Mode
        ## last value = 1 ---> VAR1 is active
        self.analyzer.write("DE CH2 ,'VX','IX',1,1")
        self.analyzer.write("DE VM1 ,'VM'")
      
        ## Turn off other channels 
        self.analyzer.write("DE CH1")
        self.analyzer.write("DE CH3")
        self.analyzer.write("DE CH4")
        self.analyzer.write("DE VS1")
        self.analyzer.write("DE VS2")
        self.analyzer.write("DE VM2")
       
    def sourceSelect(self,start, stop, step, compliance):
        ## Select [S]ource [S]elect page  
        self.analyzer.write("SS")

        ## Set up the current sweep variable. 
        ## This is a linear sweep from -0.001A to 0.001A
        CMD = "SS VR1,%s,%s,%s,%s"%(str(start), str(stop), str(step),str(compliance))
        self.analyzer.write(CMD)

        ## Set up the list display. Note that if something is not in
        ## the list we will not be able to request the data.
        self.analyzer.write("SM")
        self.analyzer.write("DM2")
        self.analyzer.write("LI 'IX','VM'")

    def measureData(self):
        ## Issue command to [M]easure [D]ata
        self.analyzer.write("MD ME1")

        ## Pause for the 4145B to populate the buffer with measurement
        ## data. If there is no sleep, Python will issue the ask commend
        ## before the measurement has been complete
        time.sleep(2)
        
        ## Print out the raw data string for inspection 
        ##print self.analyzer.ask("DO 'VSA'")
        ##print self.analyzer.ask("DO 'VSB'")
        
        ## Request data and transform it in to a list by splitting 
        ## on regular expression. Note that the slice operation is 
        ## needed because the first item in the list after re.split
        ## is the empty string.

        data0 = re.split(',*[NT]\s*', self.analyzer.ask("DO 'VM'"))[1:] 
        data1 = re.split(',*[NT]\s*', self.analyzer.ask("DO 'IX'"))[1:] 
 
        ## Cast as floats
        data0 = [float(i) for i in data0]
        data1 = [float(i) for i in data1]
      
        return [data0,data1]


if __name__ == "__main__": 

    ## A program to test things out. This can be used 
    ## to calibrate to a 50 OHM load
    GPIB = 1
    pa = hp_4145B(GPIB)
    pa.channelDefinition()
    pa.sourceSelect(-0.5,0.5,5e-2,0.01)
    [voltage,current] = pa.measureData()
    
    ## Calculate the derivitive with ndfit
    resistance =  ndfit.derivative(voltage,current)
    cal = sum(resistance)/len(resistance)
 
    fig = plt.figure(1)
    ax1 = plt.subplot(111)
    ax1.plot(voltage, current)
    ax1.set_xlabel("Voltage")
    ax1.set_ylabel("Current")
    ax1.set_ylim(-10e-3,10e-3)
    ax2 = plt.twinx(ax1) 
    h2 = ax2.plot(voltage, resistance,'r')
    ax2.set_ylabel("Resistance")
    ax2.set_ylim(40,60)
    plt.legend([h2],['R = %sOhm'%str(cal)],loc=2)

    plt.show()
