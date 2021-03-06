#!/usr/bin/env python 
import os
import sys
import visa
import ndfit
import numpy
import inspect 
import matplotlib.pyplot as plt
import cPickle as pickle

## Set the location of the data to be saved
STRUCTURE = 'q4s2s.pkl'

## Get to driver directory
os.chdir('../drivers')
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

## Actually import the driver
from OhmicDriverSeries import hp_4145B

## Move back into the current folder and set some constants 
os.chdir('../ohmic/chipA2')
size, width = [10.,20.,40.,90.],[90.]*4
squares = ndfit.quotient(size,width)

## Initialize the dude and some storage vectors
GPIB = 1
pa = hp_4145B(GPIB)
pa.channelDefinition()
data = {
    "voltage":[0]*len(size),
    "current":[0]*len(size),
    "resistance":[0]*len(size),
    "sheet":[0]*len(size),
    "size":[0]*len(size), 
    "calibration":[0]*len(size)
}


vlow, vhigh, vstep, compliance = -0.5,0.5,5e-2,0.5

## Set the current and measure the voltage
pa.sourceSelect(vlow,vhigh,vstep,compliance)
cal = 49.9048

## Start the actual measurement
i = 0
raw_input("Start Measurement: Probe Structure (enter)")
while i < len(squares):

    ## Perform the measurement
    [voltage,current] = pa.measureData()

    ## Calculate the derivitive with ndfit
    resistance =  ndfit.derivative(voltage,current)

    ## Take out the 50Ohm load
    resistance = [r-cal for r in resistance]
 
    ## set the results into the dict. Using SET rather than append
    ## avoids complications with the redo loop
    data["voltage"][i]=voltage
    data["current"][i]=current
    data["resistance"][i]=sum(resistance)/len(resistance)
    data["sheet"][i]=(sum(resistance)/len(resistance)/squares[i])
    data["size"][i]=("%sum x %sum"%(str(size[i]),str(width[i])))
    data["calibration"][i]=cal

    c = raw_input("Probe Strucute. Press (r) for redo (enter) to continue:")
    if (c == 'r'):
        i-=1
    i+=1

## Print the sheet resistance     
print data["resistance"]
print data["sheet"]
pkl_file = open(STRUCTURE,'wb')
pickle.dump(data,pkl_file)

## Plot the data
hlist = []
for i in range(len(data["size"])):
    h = plt.plot(data["voltage"][i], data["current"][i])
    hlist.append(h) 
plt.xlabel("Voltage")
plt.ylabel("Current")    
plt.legend(hlist,data["size"],loc=2)
plt.show()
