#!/usr/bin/python2.6

import math
import visa
import shlex
import time
import os

import matplotlib.pyplot as plt
# Connect to Devices

class scope_tds7104:
    def __init__(self, GPIB):
        self.scope = visa.instrument("GPIB::%s"%GPIB)

        #############################
        # Initialize Scope Channels #
        #############################
    
        # Set Horizontal Axis Parameters
        self.scope.write("HOR:SCA 1E-6")
        self.scope.write("HOR:TRIG:POS: 20%")
        # Set Horizontal Resolution 
        self.scope.write("HOR:RESO 5000")
        self.scope.write("HOR:RECO 5000")

        # Zero Vertical Offsets of Channels to Zero
        self.scope.write("CH1:OFF 0")
        self.scope.write("CH2:OFF 0")
        # Set channel copuling to DC 
        self.scope.write("CH1:COUP DC")
        self.scope.write("CH2:COUP DC")
        # Set Reasonable Vertical Scale of Channels
        self.scope.write("CH1:SCA 2")
        self.scope.write("CH2:SCA 0.02")
        self.scope.write("CH3:SCA 2")
        self.scope.write("CH4:SCA 0.02")
        
        # Turn off acquisition
        self.scope.write("ACQ:STATE OFF")
        self.scope.write("ACQ:REPE OFF")
                
        #####################
        # Set up Triggering #
        #####################

        # Trigger from AUX in. Rising Slope with DC coupling
        self.scope.write("TRIG:A:EDGE:COUP DC")
        self.scope.write("TRIG:A:EDGE:SLO RIS")
        self.scope.write("TRIG:A:EDGE:SOU AUX")


    def setScale(self, channel, scale): 
        self.scope.write("CH%s:SCA %s"%(channel, scale)) 
    
    def setSingle(self):
        # Enable Single Aquisition Mode
        self.scope.write("ACQ:STATE ON")
        self.scope.write("ACQ:REPE ON")
        self.scope.write("ACQ:STOPA SEQ")

    def setAverage(self, numavg): 
        # Enable averaging mode
        self.scope.write("ACQ:MOD AVE")
        self.scope.write("ACQ:STOPA SEQ")
        self.scope.write("ACQ:NUMAV %s"%numavg)
        self.scope.write("ACQ:STATE ON")

    def setAcqOff(self): 
        # Write to turn off aquistion this prevents 
        # averaging of sucsessive aquisitions
        self.scope.write("ACQ:STATE OFF")

    def chx_data(self, channel):
        
        # Request Encoding and curve Data
        self.scope.write("DAT:INIT")
        self.scope.write("DAT:SOU CH%s"%channel)
        self.scope.write("DAT:ENC ASCIi")
        datString = self.scope.ask("CURV?")
        
        # Obtain x and y scaling factors
        nTstep = int(self.scope.ask("HOR:RECO?"))
        tScale = float(self.scope.ask("WFMO:XIN?"))
        vScale = float(self.scope.ask("WFMO:YMU?"))

        # Build Time Data
        tList = range(0,nTstep) 
        timeData = [tScale*float(i) for i in tList] 

        # Build Votlage Data
        my_splitter = shlex.shlex(datString, posix=True)
        my_splitter.whitespace += ','
        my_splitter.whitespace_split = True
        vList = list(my_splitter)
        voltageData = [vScale*float(i) for i in vList]

        return timeData, [i-voltageData[0] for i in voltageData]

# Driver Code Here
class source_66xx:
    def __init__(self,GPIB):
        self.source = visa.instrument("GPIB::%s"%GPIB)
        self.source.ask("*IDN?")
        # Zero Voltages
        self.source.ask("VSET 1,0") 
        self.source.ask("VSET 2,0") 
        # Outputs off
        self.source.write("OUT 1,0")
        self.source.write("OUT 2,0")

    def setVoltage(self,v1set,v2set,compliance):
        if v1set > 40:
            print "V1 is too large"
        elif v2set > 40:
            print "V2 is too large"
        else:
            try:
                # Output One
                self.source.write("OUT 1,0")
                self.source.ask("VSET 1,"+str(v1set)) 
                self.source.ask("ISET 1,"+str(compliance)) 
                self.source.write("OUT 1,1")
                # Output Two
                self.source.write("OUT 2,0")
                self.source.ask("VSET 2,"+str(v2set))
                self.source.ask("ISET 2,"+str(compliance)) 
                self.source.write("OUT 2,1")
            except: 
                print "Error: Check GPIB connection and address"
    def outputsOff(self):
        self.source.write("OUT 1,0")
        self.source.write("OUT 2,0")
    def getVoltage(self):  
        print "V1="+self.source.ask("VOUT? 1")
        print "V2="+self.source.ask("VOUT? 2")
    def getCurrent(self):
        print "I1="+self.source.ask("IOUT? 1")
        print "I2="+self.source.ask("IOUT? 2")

class tabor_8500: 
 
    def __init__(self, GPIB):
        
        # Set Address 
        self.pulse = visa.instrument("GPIB::%s"%GPIB)
        # Select Channel A and turn output off
        self.pulse.write("CHA")
        self.pulse.write("D1") 
                
    def setPulseA(self,wPulse,Period):  
        # Set Parameters
        self.pulse.write("WID %sus"%wPulse) 
        self.pulse.write("PER %sms"%Period) 
        # Activate Output
        self.pulse.write("D0")
    
    def setVoltageA(self,vLow,vHigh):

        self.vLow, self.vHigh = vLow, vHigh
        if math.fabs(vHigh) > 5:
            print "Voltage is too high: output off"
            self.vHigh = 0
        
        if math.fabs(vLow) > 5:
            print "Voltage is too high: output off"
            self.vLow = 0
    
    def offPulseA(self):
        self.pulse.ask("D1")
        
    def triggerA(self):
        
        if self.vLow != self.vHigh: 
            
            self.pulse.write("LOL %sV"%round(self.vLow,2))
            self.pulse.write("HIL %sV"%round(self.vHigh,2))
            
            # Set to Triggered Mode and turn on output
            self.pulse.ask("D0")
            self.pulse.ask("M1")
            self.pulse.ask("TRG")

        else: 
            self.pulse.ask("D1")

#########
# Usage #
#########

if __name__ == "__main__":

    # Power Supply
    GPIB, vlow, vhigh, compliance = 6,0,4,.5
    source = source_66xx(GPIB)
    source.setVoltage(vlow,vhigh,compliance) 
    source.getVoltage() 
    source.getCurrent()

    # Scope
    GPIB = 15
    scope = scope_tds7104(GPIB)
    scope.setAverage(512)
    scope.setScale(1,2) 
    scope.setScale(2,0.02) 
    scope.setScale(3,2)

    # Pulse 
    GPIB = 9 
    pulse = tabor_8500(GPIB)
    pulse.setPulseA(1,1)
    pulse.setVoltageA(0,2)
    pulse.triggerA()

    # Buffer sleep
    time.sleep(1)

    # Scope Data
    tch1, vch1 = scope.chx_data(1)
    tch2, vch2 = scope.chx_data(2)
    tch3, vch3 = scope.chx_data(3)

    plt.plot(tch1,vch1)
    plt.plot(tch2,vch2)
    plt.plot(tch3,vch3)
    plt.title("Test Pulses")
    plt.xlabel("Time (s)")
    plt.ylabel("Voltage (V)")
    plt.show()
    #source.outputsOff()
