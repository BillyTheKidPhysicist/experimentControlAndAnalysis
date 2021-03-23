from tqdm import tqdm
import sys
import threading
import globalVariables as gv
import matplotlib.pyplot as plt
import time
import numpy as np
from DAQClass import DAQPin
from CameraClass import Camera
class Sweeper:
    def __init__(self,GUI):
        self.GUI=GUI
        self.cameraNear=None
        self.cameraFar=None
        self.DAQDataArr=None #array to hold data read from DAQ board. each row is a sample of data like
                # [outputVoltage, Lithium reference chamber voltage]
        self.imageArrList=[] #list to hold array of images. looks like [[imageN1,imageF1],[imageN2,imageF2],..]
            #where each image is a numpy array. So first image in each pair is for near, second for far
        self.galvoOut = DAQPin(gv.galvoOutPin)
        self.lithiumRefIn=DAQPin(gv.lithiumRefInPin)


        self.minVolt=gv.minScanVal
        self.maxVolt=gv.maxScanVal
        self.DAQVoltArr=np.linspace(self.minVolt,self.maxVolt,num=int((self.maxVolt-self.minVolt)*gv.samplesPerVoltDAQ))


        x1N=int(GUI.x1NearBox.get())
        x2N=int(GUI.x2NearBox.get())
        y1N=int(GUI.y1NearBox.get())
        y2N=int(GUI.y2NearBox.get())
        self.imageParamNear=[x1N,x2N,y1N,y2N]
        x1F=int(GUI.x1FarBox.get())
        x2F=int(GUI.x2FarBox.get())
        y1F=int(GUI.y1FarBox.get())
        y2F=int(GUI.y2FarBox.get())
        self.imageParamFar=[x1F, x2F, y1F, y2F]
        self.numExp=int(GUI.expNumBox.get())
        self.imageStartVolt=float(GUI.startVoltBox.get())
        self.imageStopVolt = float(GUI.stopVoltBox.get())
        self.expTimeNear=int(GUI.expTimeNearBox.get())
        self.expTimeFar=int(GUI.expTimeFarBox.get())
        self.binSizeNear=int(self.GUI.binSizeNearBox.get())
        self.binSizeFar=int(self.GUI.binSizeFarBox.get())



        if self.imageStartVolt<self.minVolt or self.imageStartVolt>self.maxVolt or self.imageStopVolt<self.minVolt or self.imageStopVolt>self.maxVolt:
            raise Exception('VOLTAGE RANGE FOR IMAGE SWEEP EXCEEDS MAXIMUM AND/OR MINIMUM ALLOWED RANGE!')
        if self.imageStopVolt<self.imageStartVolt:
            raise Exception('ENDING VOLTAGE IS BEFORE STARTING VOLTAGE FOR IMAGE SWEEP!')
        self.imageVoltArr=np.linspace(self.imageStartVolt,self.imageStopVolt,num=self.numExp)
    def sweep(self):
        #this sweeps the galvo output voltage. There are two arrays, DAQVoltArr and imageVoltArr. DAQVoltArr contains all
        #the voltage values to collect DAQ data at. imageVoltArr contains the values to take iamges at. imageVolt array's
        #range must be less than or equal to DAQVoltArr's range. The loop searchs for which step is next, jumps to that point
        #and then increments the counter.
        self._initialize_Cameras()
        i=0 #counter for DAQVoltArr
        j=0 #coutner for imageVoltArr
        loop=True
        volt=0
        tempList=[]
        gv.begin_Sound()
        lastImage=False #this is so that the last image is taken. It will flip from False to True once, and then no more
            #images
        takeImage=False #wether to take images
        totalSteps=self.DAQVoltArr.shape[0]+self.imageVoltArr.shape[0]

        print('\n \n \n \n')
        print('-----SWEEPING NOW----')
        time.sleep(.01) #if you don't wait a little then the progress bar and other messages will get messed up
            #in the terminal because they will try to write on top of each other
        progressBar=tqdm(total=totalSteps)
        while loop==True:
            progressBar.update()
            if i==self.DAQVoltArr.shape[0]-1 and j==self.imageVoltArr.shape[0]-1:
                loop=False
                volt = self.DAQVoltArr[i]
                if lastImage==False: #if the last image occurs at the last DAQ voltage as well
                    takeImage=True
                    lastImage=False
            else:
                if self.DAQVoltArr[i]==self.imageVoltArr[j]: #if potential next voltages are equal
                    volt = self.DAQVoltArr[i]
                    if i!=self.DAQVoltArr.shape[0]-1: #don't increment if its at the end!
                        i+=1
                    if j != self.imageVoltArr.shape[0]-1:#don't increment if its at the end!
                        j += 1
                        takeImage=True
                elif self.DAQVoltArr[i]<self.imageVoltArr[j]:
                    if i!=self.DAQVoltArr.shape[0]-1:#don't increment if its at the end!
                        volt = self.DAQVoltArr[i]
                        i+=1
                    else:
                        volt = self.imageVoltArr[i]
                        j+=1
                        takeImage=True
                elif self.imageVoltArr[j]<self.DAQVoltArr[i]:
                    if j!=self.imageVoltArr.shape[0]-1:#don't increment if its at the end!
                        volt = self.imageVoltArr[j]
                        j+=1
                        takeImage=True
                    elif lastImage==False: #special case for taking the last image
                        lastImage=True #Now it won't do this again. The loop will come here from now on, but it won't
                            #do anything but increment the galvo voltage because j!=self.imageVoltArr.shape[0]-1 will be\
                            #false and lastImage==False will be false also
                        volt = self.imageVoltArr[j]
                        takeImage=True
                    else:
                        volt = self.DAQVoltArr[i]
                        i+=1
            self.galvoOut.write(volt)
            tempList.append([volt,self.lithiumRefIn.read(numSamples=1000)])
            if takeImage==True:
                #for i in range(10):
                #    self._take_Exposures()
                #print((time.time()-t)/10)
                self.imageArrList.append(self._take_Exposures()) #the appended object is a list like [imageNear,imageFar]. If
                #there is no camera active for a given image the entry is None
                takeImage=False
        progressBar.close()
        time.sleep(.1)#like I said above. Pause to allow the progress bar to finish writting so it doesn't get messed up
        print('-----END OF SWEEP-----')
        self._close_DAQ_Pins()
        self._close_Cameras()
        self.DAQDataArr=np.asarray(tempList)
        gv.finished_Sound()
        #np.savetxt('data1.txt',self.DAQDataArr)
        #y=self.DAQDataArr[:,1]
        #plt.plot(y)
        #plt.show()
        self.GUI.imageArrList=self.imageArrList #this way if there is a previous list it is overwritten
        self.GUI.DAQDataArr=self.DAQDataArr


    def _close_DAQ_Pins(self):
        self.galvoOut.close()
        self.lithiumRefIn.close()
    def _initialize_Cameras(self):
        if self.GUI.cameraVarData.get()=='BOTH':
            self.cameraFar=Camera('FAR', self.expTimeFar, self.imageParamFar,bin=self.binSizeFar)
            self.cameraNear=Camera('NEAR', self.expTimeNear, self.imageParamNear,bin=self.binSizeNear)
        elif self.GUI.cameraVarData.get()=='NEAR':
            self.cameraNear=Camera('NEAR', self.expTimeNear, self.imageParamNear,bin=self.binSizeNear)
        elif self.GUI.cameraVarData.get()=='FAR':
            self.cameraFar=Camera('FAR',self.expTimeFar,self.imageParamFar,bin=self.binSizeFar)
        else:
            gv.error_Sound()
            raise Exception('NO VALID CAMERA NAME PROVIDED')
    def _take_Exposure_Wrapper(self,resultsDict,camera):
        #wrapper for taking images concurrently.
        if camera.camName=='NEAR':
            resultsDict['NEAR']=self.cameraNear.aquire_Image()
        elif camera.camName=='FAR':
            resultsDict['FAR']=self.cameraFar.aquire_Image()
        else:
            gv.error_Sound()
            raise Exception('NO VALID CAMERA NAME PROVIDED')
    def _take_Exposures(self):
        if self.GUI.cameraVarData.get()=='BOTH':
            resultsDict={} #this is used to add the images taken concurrently. I use a dictionary so I can keep track of
                #which image belongs to which camera
            T1=threading.Thread(target=self._take_Exposure_Wrapper,args=(resultsDict,self.cameraNear))
            T2=threading.Thread(target=self._take_Exposure_Wrapper, args=(resultsDict, self.cameraFar))
            T1.start()
            T2.start()
            T1.join()
            T2.join()
            imgNear=resultsDict['NEAR']
            imgFar=resultsDict['FAR']
            return [imgNear,imgFar]
        elif self.GUI.cameraVarData.get()=='NEAR':
            imgNear=self.cameraNear.aquire_Image()
            return [imgNear,None]
        elif self.GUI.cameraVarData.get()=='FAR':
            imgFar=self.cameraFar.aquire_Image()
            return [None,imgFar]
        else:
            gv.error_Sound()
            raise Exception('NO VALID CAMERA NAME PROVIDED')
    def _close_Cameras(self):
        if self.GUI.cameraVarData.get()=='BOTH':
            self.cameraNear.close()
            self.cameraFar.close()
        elif self.GUI.cameraVarData.get()=='NEAR':
            self.cameraNear.close()
        elif self.GUI.cameraVarData.get()=='FAR':
            self.cameraFar.close()
        else:
            raise Exception('NO VALID CAMERA NAME PROVIDED')
