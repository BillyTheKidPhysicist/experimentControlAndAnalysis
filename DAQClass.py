import nidaqmx as daq
import nidaqmx.constants as daqConstants
import nidaqmx.stream_readers as daqReader
import nidaqmx.stream_writers as daqWriter
import globalVariables as gv
import sys
import matplotlib.pyplot as plt
import time
import numpy as np
import traceback
from scipy.optimize import curve_fit


class DAQPin:
    def __init__(self,pinName,currentVolt=0.0):
        self.daqChannel=None
        self.writer=None
        self.pinName=pinName
        self.currentVolt=currentVolt #to track the current output voltage for the laser locking
        for x in range(0,len(gv.pinNameList)):
            if gv.pinNameList[x] == self.pinName:
                self.minVolt,self.maxVolt=gv.pinVoltRangeList[x]
                self.daqChannel = daq.Task()
                if gv.pinTypeList[x]=="in":
                    self.daqChannel.ai_channels.add_ai_voltage_chan("Dev1/" + self.pinName, min_val=self.minVolt, max_val=self.maxVolt)
                    self.type="in"

                if gv.pinTypeList[x]=="out":
                    self.daqChannel.ao_channels.add_ao_voltage_chan("Dev1/" + self.pinName, min_val=self.minVolt, max_val=self.maxVolt)
                    self.type="out"



        if self.daqChannel==None:
            print("------------DAQ_Error-----------------")
            time.sleep(.1), traceback.print_stack(), time.sleep(.1)
            print("you tried using a pin named \""+str(pinName)+"\"")
            print("that pin does not exist, or could not be used")
            print("see globalVariables.py for which pinNames to use")
            print("-------end DAQ_Error------------------")
            sys.exit()



    def read(self,numSamples=1,average=True,error=False):
        if numSamples==1:
            return self.daqChannel.read()
        if numSamples>1:
            reader = daqReader.AnalogSingleChannelReader(self.daqChannel.in_stream)
            dataArray = np.zeros(numSamples)
            reader.read_many_sample(dataArray, numSamples, timeout=daq.constants.WAIT_INFINITELY)
            if average==True and error==True:
                return np.average(dataArray),np.std(dataArray,ddof=1)
            if average==True and error==False:
                return np.average(dataArray)
            if average==False and error==True:
                return np.std(dataArray,ddof=1)
            else:
                return dataArray

    def write(self,value,sweep=False):
        if sweep==True or self.pinName==gv.galvoOutPin: #if writing to laser, you need to go in steps to avoid losing lock
            #self.daqChannel.write(value)
            #time.sleep(.01)
            if np.abs(value-self.currentVolt) < 1.0/gv.stepsPerVoltGalvo: #if value is very close to target, dont bother sloping
                self.daqChannel.write(value)
                self.currentVolt=value  # remembering what the output voltage is
                time.sleep(.01)
                return
            else:
                steps=int(np.abs(self.currentVolt-value)*gv.stepsPerVoltGalvo)+1 #number of steps between volts. I add one
                            #so there isn't a case of zero steps somehow.
                scanValues=np.linspace(self.currentVolt,value,steps) #the voltage value array
                waitTime=np.abs(self.currentVolt-value)*gv.timePerVolt/steps        #time to wait between steps
                for volt in scanValues:
                    self.daqChannel.write(volt)
                    time.sleep(waitTime) #wait between writing voltages
                self.currentVolt = value  # remembering what the output voltage is
                return
        else:
            self.daqChannel.write(value)
            self.currentVolt=value  # remembering what the output voltage is


    def close(self,zero=True):
        if self.type=="out":
            if zero==True:
                self.write(0.0)
        self.daqChannel.close()
        #print('DAQ board pin '+self.pinName + ' closed')


