import numpy as np
import scipy.signal as sps
import scipy.interpolate as spi
import globalVariables as gv
import matplotlib.pyplot as plt


def make_MHz_Scale(DAQData):
    temp=MHzScale(DAQData)
    temp.make_MHZ_Scale()

class MHzScale:
    def __init__(self,DAQData):
        self.DAQData=DAQData
        self.galvoVolt=self.DAQData[:, 0]  #The voltage read in from the run of the galvo output
        self.liRefVolt=self.DAQData[:, 1]  #the voltage from the diode/PMT looking at referance flourescence
    def make_MHZ_Scale(self):
        peak1VoltMean, peak2VoltMean=self.find_Peaks()
        MHzPerVolt=(gv.Li_2S_F2_F1_Sep/1E6)/(peak2VoltMean-peak1VoltMean)
        self.galvoVolt=self.DAQData[:,0]
        MHzScaleArr=(self.galvoVolt-peak1VoltMean)*MHzPerVolt
        return MHzScaleArr
    def find_Peaks(self):
        pointsPerVolt=self.galvoVolt[self.galvoVolt<self.galvoVolt[0]+1].shape[0] #about how many points per volt there are. of course
            #the camera will add more data points in the range over image aquisition occured, but that won't be at the beginning
        windowLength=int((51/102)*pointsPerVolt) # I determined this 'experimentally'
        windowLength=(windowLength//2)*2+1 #ensure window length is odd. A requirement of savgol_filter
        polyOrder=3
        liRefVoltSmooth=sps.savgol_filter(self.liRefVolt,windowLength,polyOrder) #smooth the data
        #use the smoothed data to subtract any tilt that may be present as well as offset
        y1=np.mean(liRefVoltSmooth[self.galvoVolt<self.galvoVolt[0]+.5]) #average of y values in the first 500 mv
        y2=np.mean(liRefVoltSmooth[self.galvoVolt>self.galvoVolt[-1]-.5])  #average of y values in the last 500mv
        x1=np.mean(self.galvoVolt[self.galvoVolt<self.galvoVolt[0]+.5]) #average of x values in the first 500 mv
        x2=np.mean(self.galvoVolt[self.galvoVolt>self.galvoVolt[-1]-.5]) #average of x values in the last 500 mv
        m=(y2-y1)/(x2-x1)
        yCorrection=liRefVoltSmooth[0]+m*(self.galvoVolt-self.galvoVolt[0])
        liRefVoltSmooth=liRefVoltSmooth-yCorrection # correct the data
        if np.mean(liRefVoltSmooth)<0: #typically the data is negative
            liRefVoltSmooth=-liRefVoltSmooth
        interp=spi.Rbf(self.galvoVolt,liRefVoltSmooth,smooth=1) #interpolate the data to make the data more dense. rbf is a very good
        #interpolater
        xSample=np.linspace(self.galvoVolt[0],self.galvoVolt[-1],num=10000)
        liRefVoltSmooth=interp(xSample) #now the data is more dense
        minHeight=np.max(liRefVoltSmooth)/10

        sep=int(.5*pointsPerVolt*(liRefVoltSmooth.shape[0]/self.liRefVolt.shape[0])) #seperation between peaks in points. Need
            #to scale to account for the fact that the data to be fit is now more dense
        peakPosArr, heightDict=sps.find_peaks(liRefVoltSmooth, height=minHeight, distance=sep)
        peak1VoltMean=np.mean(xSample[peakPosArr[:2]])
        peak2VoltMean=np.mean(xSample[peakPosArr[2:]])

        return peak1VoltMean, peak2VoltMean