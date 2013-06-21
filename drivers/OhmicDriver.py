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
        self.analyzer.write("DE CH2 ,'VX','IX',2,1")
        self.analyzer.write("DE VM1 ,'VS'")
 
        ## Turn off other channels 
        self.analyzer.write("DE CH1")
        self.analyzer.write("DE CH3")
        self.analyzer.write("DE CH4")
        self.analyzer.write("DE VM2")
        self.analyzer.write("DE VS1")
        self.analyzer.write("DE VS2")

    def sourceSelect(self,start, stop, step, compliance):
        ## Select [S]ource [S]elect page  
        self.analyzer.write("SS")

        ## Set up the current sweep variable. 
        ## This is a linear sweep from -0.001A to 0.001A
        CMD = "SS IR1,%s,%s,%s,%s"%(str(start), str(stop), str(step),str(compliance))
        self.analyzer.write(CMD)

        ## Pack the independent parameter into a list and return it
        current = list(np.linspace(start, stop, int((stop-start)/(step)), endpoint=False))
        current.append(stop)
        return current

    def measureData(self):
        ## Issue command to [M]easure [D]ata
        self.analyzer.write("MD ME1")

        ## Pause for the 4145B to populate the buffer with measurement
        ## data. If there is no sleep, Python will issue the ask commend
        ## before the measurement has been complete
        time.sleep(1)
        
        ## Print out the raw data string for inspection 
        print self.analyzer.ask("DO 'VS'")
        
        ## Request data and transform it in to a list by splitting 
        ## on regular expression. Note that the slice operation is 
        ## needed because the first item in the list after re.split
        ## is the empty string.
        data = re.split(',*[NT]\s*', self.analyzer.ask("DO 'VS'"))[1:] 
        return [float(i) for i in data]

if __name__ == "__main__": 

    ## A program to test things out
    GPIB = 1
    pa = hp_4145B(GPIB)
    pa.channelDefinition()
    current = pa.sourceSelect(-0.05,0.05,5e-3,7)
    voltage = pa.measureData()
    
    ## Calculate the derivitive with ndfit
    resistance =  ndfit.derivative(voltage,current)
    print sum(resistance)/len(resistance)
    plt.plot(voltage, current)
    plt.show()
