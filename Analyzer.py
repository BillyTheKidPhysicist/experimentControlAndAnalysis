import time
import sys
import numpy as np
import scipy.optimize as spo
import scipy.special as spf
import scipy.signal as spsi
import globalVariables as gv

import matplotlib.pyplot as plt

class Analyzer:
    def __init__(self,imagesSumArr,imageFreqMHzArr,S):
        #imagesSumArr: array where each entry is the sum of the pixels in that image over a region
        #imageFreqMHzArr: MHz values that correspond to images in imagesSumArr
        #S: I/I_sat. saturation constant
        self.imagesSumArr=imagesSumArr
        self.imageFreqMHzArr=imageFreqMHzArr
        self.S=S
        self.fitFunc=None #A function that returns the spectral profile of the fit
        self.F0=None #the value of the center of the peak of the profile
        self.FWHM=None #the FWHM of the data
        self.T=None #the temperature of the fit




    def fit_Image_Data(self):
        #   Billy Debenham,Jly 3rd 2018. billydebenham@gmail.com
        #   DESCRIPTION:
        #   This method returns the fitting parameters of the triple voigt
        #
        #   how it works is the following
        #   1.Guess paramters like the height,y axis offset and x axis offset
        #   2.fit the data
        #
        #   INPUTS:
        #   MHzScale: a 1D array MHZ relative to the first lambdip peak
        #   pixelValue: a 1D array of intensity of a pixel or pixels
        #
        #
        #   NOTES:
        #   If you want to change the guess parameter and bounds of the other variables,
        #   find the 3 arrays that hold that information and change it there.
        #   MHzScale and pixelValue must have same dimension
        #
        #   KEY VARIABLES
        #   a: height scaling value
        #   b: center value of frequency
        #   c: offset from zero
        #   dF: FWHM of doppler broadening
        #
        #   RETURNS
        #   PF: array of paramters in voigt fit
        #   stdDev: array of standard deviations of parameters



        #func=spi.interp1d(MHzScale,pixelValue,kind='quadratic')
        #func=spi.Rbf(MHzScale,pixelValue)
        #xNew=np.linspace(MHzScale[0],MHzScale[-1],10000)
        #pixelValueSD=func(xNew)

        #plt.plot(MHzScale,pixelValue)
        #plt.plot(xNew,func(xNew))
        #plt.show()


        #smooth=self.imagesSumArr#spsi.savgol_filter(pixelValue,pixelValue.shape[0]//10+1,3) temporarily changed
        #minHeight=(np.max(smooth)-np.min(smooth))/2.0 +np.min(smooth)
        #peak=spsi.find_peaks(smooth,height=minHeight)[0][0]
        error=np.ones(self.imageFreqMHzArr.shape[0])
        #error[peak+2]=.2
        #error[peak-2]=.2
        #error[peak+1]=.1
        #error[peak-1]=.1
        #error[peak]=.01

        #--------GUESSING-----------------------
        #aG: height constant
        #bG: center frequency
        #cG: vertical offset
        #dFG: FWHM of doppler broadening
        aG=np.max(self.imagesSumArr)-np.min(self.imagesSumArr)    #height constant
        bG=self.imageFreqMHzArr[np.argmax(self.imagesSumArr)]  #x axis offset
        cG=np.average(self.imagesSumArr[-int(self.imagesSumArr.shape[0]/10.0):]) #Set the y axisoffsest guess based on the 
                    # average of the last 10% of data
        dFG=20  #good guess for FWHM of gaussian

        y1=np.mean(self.imagesSumArr[:5]) #take average of first few data points
        y2=np.mean(self.imagesSumArr[-5:]) #take the average of the last few data points
        x1=np.mean(self.imageFreqMHzArr[:5])
        x2=np.mean(self.imageFreqMHzArr[-5:])
        tiltmG=(y2-y1)/(x2-x1)


        guess = np.array([aG    ,bG            , cG*1.0, dFG ,tiltmG]) #array of guess

        PF, pcov = spo.curve_fit(self.spectral_Profile_Wrapper, self.imageFreqMHzArr, self.imagesSumArr, p0=guess)
        time.sleep(.1)
        perr=np.sqrt(np.diag(pcov))
        print(PF)
        print(perr)



        self.F0=self.find_Voigt_Peak(self.imageFreqMHzArr,PF)
        dF=PF[3] #FWHM of spectral profile

        def temp(freq):
            #freq: HHz
            return self.spectral_Profile_Wrapper(freq, *PF)
        self.fitFunc=temp
        print(dF)
        self.T=self.calculateTemp(dF)
        print(self.T*1E3)
        plt.scatter(self.imageFreqMHzArr, self.imagesSumArr, c='red',marker='x')
        x=np.linspace(self.imageFreqMHzArr[0],self.imageFreqMHzArr[-1],num=1000)
        plt.plot(x, self.spectral_Profile_Wrapper(x, *PF))
        plt.show()
        


#---------SINGLE VOIGT FIT---------
#this voigt is normalized to a height of 1, then multiplied by the variable a
#there is no yaxis offset in the voigt, that is done later
    def voigt(self,f,a,f0,dF,gamma):
        #units must be consistent!!
        #f, frequency value
        #a, height of voigt
        #f0, center frequency
        #dF, FWHM of the gaussian
        #gamma, FWHM of the lorentzian
        sigma=dF/(2*np.sqrt(2*np.log(2))) #convert gaussian FWHM to std
        gamma=gamma/2 #convert lorentzian FWHM to HWHM
        x=f-f0
        v0=spf.voigt_profile(0,sigma,gamma)
        v=a*spf.voigt_profile(x,sigma,gamma)/v0
        return v
    def spectral_Profile_Wrapper(self,freq,a,b,c,dF,tiltm):
        #this is wrapper for the sextuple voigt fit that also adds an offset, and maybe in the future the ability to
        #compensate for saturation intensity
        gamma=gv.LiTau/1E6
        val= self.spectral_Profile(freq, a, b, c, dF, gamma)
        val=val+(tiltm*(freq-self.imageFreqMHzArr[0]))
        return val
#------SEXTUPLE VOIGT FIT----------
    def spectral_Profile(self, freq, a, b, c, dF,gamma ):
        #this make a sextuple voigt fit
        #units must be consitant!
        #freq: frequency
        #a: height constant
        #b: center frequency
        #c: vertical offset
        #dF: FWHM of doppler broadening
        #gamma: FWHM of lorentzian

        #gamma=gamma*np.sqrt(1+S) #linewidth broadening

        ratio=4*2.54#ratio of intensity of f=2 to f=1. First number is power ratio in bands. Second fraction is ratio
                    #of peaks in lithium reference cell
        a2=ratio*a
        a1=a
        val=0
        ##F=2 transition
        val+=a2*self.voigt(freq, gv.S21, b+gv.F1Sep/1E6, dF, gamma)
        val+=a2*self.voigt(freq, gv.S22, b+gv.F2Sep/1E6, dF, gamma)
        val+=a2*self.voigt(freq, gv.S23, b+gv.F3Sep/1E6, dF, gamma)
#
#
        ##F=1 transition
        val+=a1*self.voigt(freq, gv.S10, b+gv.F0Sep/1E6, dF, gamma)
        val+=a1*self.voigt(freq, gv.S11, b+gv.F1Sep/1E6, dF, gamma)
        val+=a1*self.voigt(freq, gv.S12, b+gv.F2Sep/1E6, dF, gamma)
        val+=c #add the offset
        return val

    def calculateTemp(self,dF,f0=gv.Li_D2_Freq,MHz=True):
        #calculate the temperatute in mk for a given doppler width at a given frequency F0
        if MHz==True:
            dF=dF*1E6
        return (dF/f0)**2 *(gv.massLi7*gv.cLight**2)/(8*gv.kB*np.log(2))

    def cheatData(self,offset=5,showFit=True,returnFit=True):

        #need to include better interaction with the plot! doesn't work as expected when accessed from mainGUI


        smooth=spsi.savgol_filter(self.imagesSumArr,21,5)
        minHeight=(np.max(smooth)-np.min(smooth))/2.0 +np.min(smooth)
        peak=spsi.find_peaks(smooth,height=minHeight)[0][0]

        x=np.flip(self.imageFreqMHzArr[peak-offset:])
        y=np.flip(self.imagesSumArr[peak-offset:])
        #newPeak=184
        error=np.ones(x.shape[0])

        PF,stdDev=self.fit_Image_Data(x,y)
        a, b, c, dF = PF
        if showFit==True:
            voigtFit=[]
            for elem in x:
                voigtFit.append(self.spectral_Profile(elem, a, b, c, dF))

            plt.title("\'Cheat\' method. MAKE SURE FIT LOOKS GOOD")
            plt.xlabel("pixel value")
            plt.ylabel("MHz")
            plt.plot(x,y)
            plt.plot(x,voigtFit)
            plt.ion()
            plt.show()
            plt.pause(.001)
            input('if the above graph looks good, click enter to continue')
            plt.ioff()
            plt.close('all')
            print("left")

        if returnFit==True:
            voigtFit=np.asarray(self.spectral_Profile(self.imageFreqMHzArr, a, b, c, dF))
            return [a, b, c, dF], stdDev, voigtFit
        return PF, stdDev

    def find_Voigt_Peak(self,imageFreqMHzArr,PF):
        x=np.linspace(self.imageFreqMHzArr[0],self.imageFreqMHzArr[-1],num=1000000)
        y=self.spectral_Profile_Wrapper(x, *PF)
        F0=x[spsi.find_peaks(y)[0]]
        return F0[0]






