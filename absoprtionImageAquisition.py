import sys
import time
from CameraClass import Camera
import numpy as np
from DAQClass import DAQPin
import globalVariables as gv
from astropy.io import fits
import os



class AbsorptionImager:
    def __init__(self,path,runName,voltStart,voltEnd,numPoints,exposure_Seconds):
        self.path=path
        self.runName=runName
        self.voltStart=voltStart
        self.numPoints=numPoints
        self.voltArr=np.linspace(voltStart, voltEnd, numPoints)
        self.exposure_MilliSeconds=exposure_Seconds*1000

        self.MKSFullScale=500.0  #full scale of N2, sccm
        self.scaleFact=1.4  #scael factor for He
        self.flowRate=50.0
        self.binning=1 #don't use binnig for absorption



        self.nozzleWaitTime_Seconds=5 #time to wait for turning on or off the nozzle
        self.waitTimeToFlow_Seconds=self.nozzleWaitTime_Seconds-exposure_Seconds
        if self.waitTimeToFlow_Seconds<0.0:
            self.waitTimeToFlow_Seconds=0.0
        self.waitTimeToStopFlow_Seconds=self.nozzleWaitTime_Seconds

        self.camera=Camera('FAR',self.exposure_MilliSeconds,bin=self.binning)
        self.galvoOut=DAQPin(gv.galvoOutPin)
        self.shutterOut=DAQPin(gv.shutterPin) #open the shutter contorl pin
        self.darkImageList=[]
        self.noFlowImageList=[]
        self.flowImageList=[]
        self.absorptionSignalImageList=[]


    def catch_Errors(self):
        if self.flowRate<=0:
            gv.error_Sound()
            raise Exception('The flowrate is zero or invalid!')
        if self.exposure_MilliSeconds<500:
            gv.error_Sound()
            raise Exception('The camera exposure is too low. Remember you enter it as seconds')
    def make_Info_File(self):
        with open('info.txt','w') as file:
            file.write('Run Info \n')
            file.write('Exposure time: '+str(self.exposure_MilliSeconds)+' ms \n')
            file.write('flow wait time: '+str(self.nozzleWaitTime_Seconds)+'s \n')
            file.write('binning :'+str(self.binning)+'\n')
            file.write('Nozzle flow rate: '+str(self.flowRate)+' SCCM \n')
            file.write('Image galvo voltages: '+str(self.voltArr)+' volts \n')
    def save_Images(self,name,imageList):
        hdu=fits.PrimaryHDU(np.asarray(imageList))
        hdul=fits.HDUList([hdu])
        hdul.writeto(name+'.fits')
    def make_Flow(self,flowDesired):
        flowOut=DAQPin(gv.flowOutPin)
        if flowDesired>500.0:
            raise Exception('REQUESTED FLOW IS GREATER THAN MAXIMUM')
        elif flowDesired>0.0:
            volt=(flowDesired/(self.MKSFullScale*self.scaleFact))*5.0
            flowOut.write(volt)
        else:
            flowOut.write(0.0)
        flowOut.close(zero=False)

    def take_Dark_Image(self):
        self.shutterOut.close_Shutter()
        image=self.camera.aquire_Image().astype(float)
        self.shutterOut.open_Shutter()
        return image
    def take_No_Flow_Image(self):
        image=self.camera.aquire_Image().astype(float)
        return image
    def take_Flow_Image(self):
        image=self.camera.aquire_Image().astype(float)
        return image
    def make_Copies_To_Not_Modify_Originals(self,darkImage, noFlowImage, flowImage):
        #copy the images to prevent modifying them by accident!
        noFlowImage=noFlowImage.copy()
        flowImage=flowImage.copy()
        darkImage=darkImage.copy()
        return darkImage,noFlowImage,flowImage
    def construct_Absorption_Signal_Image(self,darkImage, noFlowImage, flowImage):
        eps=1.0 #to avoid divide by zero
        darkImage,noFlowImage,flowImage=self.make_Copies_To_Not_Modify_Originals(darkImage,noFlowImage,flowImage)
        noFlowImage=noFlowImage-darkImage
        flowImage=flowImage-darkImage
        flowImage[np.abs(flowImage)<eps]=eps #get rid of small numbers to prevent divide by zero
        noFlowImage[np.abs(noFlowImage)<eps]=eps
        signalImage=flowImage/noFlowImage
        return signalImage
    def save_Image_Lists(self):
        self.save_Images(self.runName+'_darkImages',self.darkImageList)
        self.save_Images(self.runName+'_noFlowImages',self.noFlowImageList)
        self.save_Images(self.runName+'_flowImages',self.flowImageList)
        self.save_Images(self.runName+'_absorptionSignalImages',self.absorptionSignalImageList)
    def change_Directory_And_Catch_File_Errors(self):
        imagesFolder=self.runName+'Folder'
        try:
            os.mkdir(self.path+'\\'+imagesFolder)
        except FileExistsError:
            gv.error_Sound()
            print('That file already exists!!')
            exit()
        except: raise Exception('some other file issue')
        os.chdir(self.path+'\\'+imagesFolder)
    def update_Image_Lists(self,darkImage,noFlowImage,flowImage,absorptionSignalImage):
        self.darkImageList.append(darkImage)
        self.noFlowImageList.append(noFlowImage)
        self.flowImageList.append(flowImage)
        self.absorptionSignalImageList.append(absorptionSignalImage)
    def close_Pins(self):
        self.galvoOut.close()
        self.shutterOut.close()
    def wait_To_Flow_After_Dark_Image(self,):
        time.sleep(self.waitTimeToFlow_Seconds)
    def wait_To_Stop_Flow_After_Flow_Image(self):
        time.sleep(self.waitTimeToStopFlow_Seconds)
    def run(self):
        self.change_Directory_And_Catch_File_Errors()
        self.catch_Errors()
        gv.begin_Sound()
        self.make_Flow(0.0)
        self.shutterOut.open_Shutter()
        for volt in self.voltArr:
            print(volt)
            self.galvoOut.write(volt)
            noFlowImage=self.take_No_Flow_Image()
            self.make_Flow(self.flowRate)
            darkImage=self.take_Dark_Image()
            self.wait_To_Flow_After_Dark_Image()
            flowImage=self.take_Flow_Image()
            self.make_Flow(0.0)
            self.wait_To_Stop_Flow_After_Flow_Image()
            absorptionSignalImage=self.construct_Absorption_Signal_Image(darkImage, noFlowImage, flowImage)
            self.update_Image_Lists(darkImage,noFlowImage,flowImage,absorptionSignalImage)
        self.close_Pins()
        self.save_Image_Lists()
        self.make_Info_File()
        gv.finished_Sound()

folderPath='C:\Data\Runs\9_23_21'
runName='run27'
voltStart=-.80
voltEnd=-.78
numExposure=4
exposureTime_Seconds=2.0
absorptionImager=AbsorptionImager(folderPath, runName, voltStart, voltEnd, numExposure, exposureTime_Seconds)
absorptionImager.run()