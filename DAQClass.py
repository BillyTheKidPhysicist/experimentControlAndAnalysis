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
import threading
from scipy.optimize import curve_fit


class DAQPin:
    def __init__(self,pinName,currentVolt=0.0):
        self.daqChannel=None
        self.channelObject =None #the object created when adding a daq channel
        self.writer=None
        self.pinName=pinName
        self.thread=None
        self.currentVolt=currentVolt #to track the current output voltage for the laser locking
        for x in range(0,len(gv.pinNameList)):
            if gv.pinNameList[x] == self.pinName:
                self.minVolt,self.maxVolt=gv.pinVoltRangeList[x]
                self.daqChannel = daq.Task()
                if gv.pinTypeList[x]=="in":
                    self.channelObject=self.daqChannel.ai_channels.add_ai_voltage_chan("Dev1/" + self.pinName,
                                                                min_val=self.minVolt, max_val=self.maxVolt)
                    self.type="in"
                if gv.pinTypeList[x]=="out":
                    self.channelObject=self.daqChannel.ao_channels.add_ao_voltage_chan("Dev1/" + self.pinName,
                                                                min_val=self.minVolt, max_val=self.maxVolt)
                    self.type="out"
                if gv.pinTypeList[x]=='counterOut':
                    self.channelObject=self.daqChannel.co_channels.add_co_pulse_chan_freq("Dev1/"+gv.servoPin) #sets to
                    # default frequency and duty cycle
                    self.daqChannel.timing.cfg_implicit_timing(
                        sample_mode=daqConstants.AcquisitionType.CONTINUOUS)
                    self.type="counterOut"
                if gv.pinTypeList[x]=='digitalOut':
                    self.channelObject=self.daqChannel.do_channels.add_do_chan("Dev1/"+gv.shutterPin)
                    self.type="digitalOut"


        if self.daqChannel==None:
            print("------------DAQ_Error-----------------")
            time.sleep(.1), traceback.print_stack(), time.sleep(.1)
            print("you tried using a pin named \""+str(pinName)+"\"")
            print("that pin does not exist, or could not be used")
            print("see globalVariables.py for which pinNames to use")
            print("-------end DAQ_Error------------------")
            sys.exit()


    def generate_Digital_Pulse_Train(self,freq=None,duty=None,width=None,sep=None):
        #provided units must be SI
        if self.pinName!=gv.servoPin:
            print("-------------ERROR----------------------------")
            print('YOU ARE TRYING TO SEND A DIGITAL PULSE THROUGH THE WRONG PIN')
            print("----------------end ERROR-------------")
            gv.error_Sound()
            sys.exit()
        #either freq or duty must not be None, or width and sep must not be None
        error=False
        if freq is not None and width is not None:
            error=True
        if duty is not None and sep is not None:
            error=True
        if freq is None and sep is None:
            error=True
        if duty is None and width is None:
            error=True
        if error==True:
            print("-------------ERROR----------------------------")
            print('YOU CAN ONLY SET PULSE WIDTH AND SEPERATION OR FREQUENCY AND DUTY')
            print("----------------end ERROR-------------")
            gv.error_Sound()
            sys.exit()
        if sep is not None and width is not None:
            if (width<sep)==False:
                print("-------------ERROR----------------------------")
                print('PULSE WIDTH MUST BE LESS THAN SEPERATION')
                print("----------------end ERROR-------------")
                gv.error_Sound()
                sys.exit()
            freq=1/sep
            duty=width/sep
        if self.thread is not None:
            #if the channel and thread is already running it needs to be stopped
            self.thread.join()
            self.daqChannel.stop()
        self.channelObject.co_pulse_duty_cyc=duty
        self.channelObject.co_pulse_freq=freq
        self.thread=threading.Thread(target=self.daqChannel.start())
        self.thread.start()

    def read(self,numSamples=1,average=True,error=False):
        #numsamples: number of samples for DAQboard to read in
        #average: wether to return the average value of samples read in rather than the samples
        #error: Return the standard error of the samples instead of the samples
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
    def write_High(self):
        #write a logical high value, ie 5 volts
        self.daqChannel.write(True)
    def write_Low(self):
        #write a logical low value, ie 0 volts
        self.daqChannel.write(False)
    def write(self,value,sweep=False):
        #value: voltage value to write to the galvo
        if sweep==True or self.pinName==gv.galvoOutPin: #if writing to laser, you need to go in steps to avoid losing lock
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

    def set_Servo_Position(self,deg):
        if self.pinName!=gv.servoPin:
            print("-------------ERROR----------------------------")
            print('YOU CAN ONLY CONTROL THE SERVO THROUGH THE SERVO PIN')
            print("----------------end ERROR-------------")
            gv.error_Sound()
            sys.exit()
        m=.5/90 #ms/degree
        b=1.5 #ms
        width= m*deg+b
        width=width/1e3 #convert to seconds
        sep=20E-3 #pulse seperation, seconds    az CVBNM,.CVBNM
        self.generate_Digital_Pulse_Train(width=width,sep=sep)



    def close(self,zero=True):
        if self.thread is not None:
            self.thread.join()
        if self.type=="out":
            if zero==True:
                self.write(0.0)
        if self.type=="digitalOut":
            if zero==True:
                self.write_Low()
        self.daqChannel.stop()
        self.daqChannel.close()
