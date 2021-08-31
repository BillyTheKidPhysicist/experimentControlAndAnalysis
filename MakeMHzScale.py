import numpy as np
import scipy.signal as sps
import scipy.interpolate as spi
import globalVariables as gv
import matplotlib.pyplot as plt
import scipy.optimize as spo
import sys
import scipy.special as spspec
import numba

def make_MHz_Scale(DAQData,returnFitFunc=False):
    MHZScaleMaker=MHzScale(DAQData)
    scale=MHZScaleMaker.make_MHZ_Scale(returnFitFunc)
    return scale
@numba.njit()
def gaussian_Normalized_Height(x,x0,sigma):
    return np.exp(-.5*((x-x0)/sigma)**2)

class MHzScale:
    def __init__(self,DAQData):
        self.DAQData=DAQData
        self.galvoVoltArrRaw=self.DAQData[:, 0].copy()  #The voltage read in from the run of the galvo output
        self.liRefVoltArrRaw=self.DAQData[:, 1].copy()  #the voltage from the diode/PMT looking at referance flourescence
        self.LiRefVoltArr=None
    def make_MHZ_Scale(self,returnFitFunc):
        #assumes the frequency versus voltage curve is straight
        peak1VoltMean, peak2VoltMean,fitFuncVolt=self.get_Peaks_And_Fit_Function()
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
        self.galvoVoltArrRaw=self.DAQData[:,0]
        MHzScaleArr=(self.galvoVoltArrRaw-peak1VoltMean)*MHzPerVolt
        # plt.plot(MHzScaleArr,fitFuncVolt(self.galvoVoltArrRaw))
        # plt.show()
        fitFunc= lambda x: fitFuncVolt((x/MHzPerVolt+peak1VoltMean)) #convert voltage to frequency for fit func
        if returnFitFunc==False:
            return MHzScaleArr
        else:
            return MHzScaleArr,fitFunc
    @staticmethod
    @numba.njit()
    def signal_Profile(v,v0, doubletSepMajor,doubletSepMinor,  aPeak1, aPeak2, aPeak3, aPeak4, aLeft,aCenter,aRight,
                       sigmaMain,sigmaLeft,sigmaCenter,sigmaRight,offset):

        signal=aPeak1*gaussian_Normalized_Height(v,v0-doubletSepMinor/2,sigmaMain)
        signal+=aPeak2*gaussian_Normalized_Height(v,v0+doubletSepMinor/2,sigmaMain)
        signal+=aPeak3*gaussian_Normalized_Height(v,v0+doubletSepMajor-doubletSepMinor/2,sigmaMain)
        signal+=aPeak4*gaussian_Normalized_Height(v,v0+doubletSepMajor+doubletSepMinor/2,sigmaMain)
        signal+=aCenter*gaussian_Normalized_Height(v,v0+doubletSepMajor/2,sigmaCenter)
        signal+=aLeft*gaussian_Normalized_Height(v,v0,sigmaLeft)
        signal+=aRight*gaussian_Normalized_Height(v,v0+doubletSepMajor,sigmaRight)
        signal+=offset
        return signal
    def get_Model_Bounds(self):
        a1And2Min=np.max(self.LiRefVoltArr)*.8
        a1And2Max=np.max(self.LiRefVoltArr)*1.2
        a3And4Min=np.max(self.LiRefVoltArr)*.2
        a3And4Max=1.0*np.max(self.LiRefVoltArr)

        aRightMin=-1.0*np.max(self.LiRefVoltArr)
        aRightMax=0.0

        aCenterMin=0.0
        aCenterMax=1.0*np.max(self.LiRefVoltArr)

        aLeftMin=-1.0*np.max(self.LiRefVoltArr)
        aLeftMax=0.0

        sigmaPeakMin=0.0
        sigmaPeakMax=.5

        sigmaLeftDipMin=0.05
        sigmaLeftDipMax=.5

        sigmaCenterMin=0.05
        sigmaCenterMax=.5

        sigmaRightDipMin=0.05
        sigmaRightDipMax=.5

        F0Min=self.galvoVoltArrRaw.min()
        F0Max=self.galvoVoltArrRaw.max()

        doubletMinorSepMin=.35
        doubletMinorSepMax=.7

        doubletMajorSepMin=1.2
        doubletMajorSepMax=1.6

        offsetMin=-.1
        offsetMax=.1


        bounds=[(F0Min,F0Max),(doubletMajorSepMin,doubletMajorSepMax),(doubletMinorSepMin,doubletMinorSepMax)]
        bounds.extend([(a1And2Min,a1And2Max),(a1And2Min,a1And2Max),(a3And4Min,a3And4Max),(a3And4Min,a3And4Max)])
        bounds.extend([(aLeftMin,aLeftMax),(aCenterMin,aCenterMax),(aRightMin,aRightMax)])
        bounds.extend([(sigmaPeakMin,sigmaPeakMax),(sigmaLeftDipMin,sigmaLeftDipMax),(sigmaCenterMin,sigmaCenterMax)])
        bounds.extend([(sigmaRightDipMin,sigmaRightDipMax),(offsetMin,offsetMax)])
        return bounds

    def remove_Tilt(self,signalArr):
        voltageWidthForMean=1.0
        y1=np.mean(signalArr[self.galvoVoltArrRaw<self.galvoVoltArrRaw[0]+voltageWidthForMean])
        y2=np.mean(signalArr[self.galvoVoltArrRaw>self.galvoVoltArrRaw[-1]-voltageWidthForMean])
        x1=np.mean(self.galvoVoltArrRaw[self.galvoVoltArrRaw<self.galvoVoltArrRaw[0]+voltageWidthForMean])
        x2=np.mean(self.galvoVoltArrRaw[self.galvoVoltArrRaw>self.galvoVoltArrRaw[-1]-voltageWidthForMean])
        m=(y2-y1)/(x2-x1)
        yCorrection=signalArr[0]+m*(self.galvoVoltArrRaw-self.galvoVoltArrRaw[0])
        signalArrFlat=signalArr-yCorrection  # correct the data
        return signalArrFlat
    def model_Fit_Cost(self,args):
        cost=np.sqrt(np.sum((self.LiRefVoltArr-self.signal_Profile(self.galvoVoltArrRaw, *args))**2))
        return cost
    
    def get_Peaks_And_Fit_Function(self):
        self.LiRefVoltArr=self.liRefVoltArrRaw.copy()
        self.LiRefVoltArr=-self.LiRefVoltArr #flip the data, it is 'upside down' from the PMT
        self.LiRefVoltArr-=self.LiRefVoltArr.min()
        self.LiRefVoltArr=self.remove_Tilt(self.LiRefVoltArr)
        self.LiRefVoltArr=self.LiRefVoltArr/self.LiRefVoltArr.max()

        import time


        bounds=self.get_Model_Bounds()

        np.random.seed(42) #to make results repeatable
        sol=spo.differential_evolution(self.model_Fit_Cost,bounds,mutation=.5)
        np.random.seed(int(time.time()))
        params=sol.x
        modelFitCostMaximum=.001
        if sol.fun/self.LiRefVoltArr.shape[0]>modelFitCostMaximum:
            plt.plot(self.galvoVoltArrRaw, self.LiRefVoltArr)
            plt.plot(self.galvoVoltArrRaw, self.signal_Profile(self.galvoVoltArrRaw, *params))
            plt.show()
            raise Exception('Model appears to fit poorly. This may overly convservative though')
        params=sol.x
        fitFunc=lambda x: self.signal_Profile(x,*params)

        peak1Volt=params[0]
        peak2Volt=peak1Volt+params[1]
        # plt.plot(self.galvoVoltArrRaw, self.LiRefVoltArr)
        # plt.plot(self.galvoVoltArrRaw, self.signal_Profile(self.galvoVoltArrRaw, *params))
        # plt.axvline(x=peak1Volt,c='black')
        # plt.axvline(x=peak2Volt,c='black')
        # plt.show()
        return peak1Volt,peak2Volt,fitFunc