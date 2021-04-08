import numpy as np
import scipy.signal as sps
import scipy.interpolate as spi
import globalVariables as gv
import matplotlib.pyplot as plt


def make_MHz_Scale(DAQData,returnFitFunc=False):
    MHZScaleMaker=MHzScale(DAQData)
    scale=MHZScaleMaker.make_MHZ_Scale(returnFitFunc)
    return scale

class MHzScale:
    def __init__(self,DAQData):
        self.DAQData=DAQData
        self.galvoVolt=self.DAQData[:, 0]  #The voltage read in from the run of the galvo output
        self.liRefVolt=self.DAQData[:, 1]  #the voltage from the diode/PMT looking at referance flourescence
    def make_MHZ_Scale(self,returnFitFunc):
        #assumes the frequency versus voltage curve is straight
        peak1VoltMean, peak2VoltMean,fitFuncVolt=self.find_Peaks()
        MHzPerVolt=(gv.Li_2S_F2_F1_Sep/1E6)/(peak2VoltMean-peak1VoltMean)
        if not gv.typicalMHzPerVolt*.95<MHzPerVolt<1.05*gv.typicalMHzPerVolt: #if the scale seems suspicous
            gv.warning_Sound(noWait=True)
            print('------------WARNING----------------')
            print("The MHz per volt is out of the expected range.")
            print("Typical value is "
                          +str(int(gv.typicalMHzPerVolt))
                          +'. \nCurrent value is '+str(int(MHzPerVolt))+'.')
            print('Its possible nothing is wrong')
            print('It is not known what the \'expected\' range is at this moment')
            print('---------END WARNING---------------')
        self.galvoVolt=self.DAQData[:,0]
        MHzScaleArr=(self.galvoVolt-peak1VoltMean)*MHzPerVolt
        fitFunc= lambda x: fitFuncVolt((x/MHzPerVolt+peak1VoltMean)) #convert voltage to frequency for fit func
        # plt.plot(MHzScaleArr,fitFunc(MHzScaleArr))
        # plt.axvline(x=0)
        # plt.grid()
        # plt.show()
        if returnFitFunc==False:
            return MHzScaleArr
        else:
            return MHzScaleArr,fitFunc
    def find_Peaks(self):
        pointsPerVolt=self.galvoVolt[self.galvoVolt<self.galvoVolt[0]+1].shape[0] #about how many points per volt there are. of course
            #the camera will add more data points in the range over image aquisition occured, but that won't be at the beginning
        windowLength=int((21/102)*pointsPerVolt) # I determined this 'experimentally'. This depends on the width
        #of the peaks
        windowLength=(windowLength//2)*2+1 #ensure window length is odd. A requirement of savgol_filter
        polyOrder=3
        liRefVoltSmooth=sps.savgol_filter(self.liRefVolt,windowLength,polyOrder) #smooth the data
        #use the smoothed data to subtract any tilt that may be present as well as offset\
        y1=np.mean(liRefVoltSmooth[self.galvoVolt<self.galvoVolt[0]+.5]) #average of y values in the first 500 mv
        y2=np.mean(liRefVoltSmooth[self.galvoVolt>self.galvoVolt[-1]-.5])  #average of y values in the last 500mv
        x1=np.mean(self.galvoVolt[self.galvoVolt<self.galvoVolt[0]+.5]) #average of x values in the first 500 mv
        x2=np.mean(self.galvoVolt[self.galvoVolt>self.galvoVolt[-1]-.5]) #average of x values in the last 500 mv
        m=(y2-y1)/(x2-x1)
        yCorrection=liRefVoltSmooth[0]+m*(self.galvoVolt-self.galvoVolt[0])
        liRefVoltSmooth=liRefVoltSmooth-yCorrection # correct the data
        if np.mean(liRefVoltSmooth)<0: #typically the data is negative
            liRefVoltSmooth=-liRefVoltSmooth
        fitFunc=spi.Rbf(self.galvoVolt,liRefVoltSmooth,smooth=1.0) #interpolate the data to make the data more dense. rbf is a very good
        #interpolater
        xSample=np.linspace(self.galvoVolt[0],self.galvoVolt[-1],num=10000)
        liRefVoltSmoothDense=fitFunc(xSample) #now the data is more dense
        minHeight=np.max(liRefVoltSmoothDense)/10

        sep=int(.5*pointsPerVolt*(liRefVoltSmoothDense.shape[0]/self.liRefVolt.shape[0])) #seperation between peaks in points. Need
            #to scale to account for the fact that the data to be fit is now more dense
        peakPosArr, heightDict=sps.find_peaks(liRefVoltSmoothDense, height=minHeight, distance=sep)
        peak1VoltMean=np.mean(xSample[peakPosArr[:2]])
        peak2VoltMean=np.mean(xSample[peakPosArr[2:]])
        # plt.plot(self.galvoVolt,-(self.liRefVolt-yCorrection),label='data')
        # plt.plot(xSample,liRefVoltSmoothDense,label='fit')
        # plt.axvline(x=peak1VoltMean)
        # plt.axvline(x=peak2VoltMean)
        # plt.legend()
        # plt.show()

        return peak1VoltMean, peak2VoltMean,fitFunc