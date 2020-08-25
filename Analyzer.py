import numpy as np
import scipy.optimize as spo
import scipy.special as spf
import scipy.signal as spsi
import globalVariables as gv

import matplotlib.pyplot as plt

class Analyzer:
    def __init__(self,imagesAvgArr,MHzScaleArr):
        self.imageAvgArr=imagesAvgArr
        self.MHzScaleArr=MHzScaleArr
        self.fitFunc=None #A function that returns the spectral profile of the fit
        self.F0=None #the value of the center of the peak of the profile
        self.FWHM=None #the FWHM of the data
        self.T=None #the temperature of the fit


    def gauss(self,T,v,m=gv.massLi7,v0=0):
        t1=np.sqrt(m/(2*np.pi*gv.kB*T))
        t2=np.exp(-m*np.power(v-v0,2)/(2*gv.kB*T))
        return t1*t2


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
        #   dF0: FWHM of doppler broadening
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


        smooth=self.imageAvgArr#spsi.savgol_filter(pixelValue,pixelValue.shape[0]//10+1,3) temporarily changed
        minHeight=(np.max(smooth)-np.min(smooth))/2.0 +np.min(smooth)
        peak=spsi.find_peaks(smooth,height=minHeight)[0][0]
        error=np.ones(self.MHzScaleArr.shape[0])
        #error[peak+2]=.2
        #error[peak-2]=.2
        #error[peak+1]=.1
        #error[peak-1]=.1
        #error[peak]=.01

        #--------GUESSING-----------------------
        #aG: height constant
        #bG: center frequency
        #cG: vertical offset
        #dF0G: FWHM of doppler broadening
        aG=np.max(self.imageAvgArr[10:-10])-np.average(self.imageAvgArr[10:-10])    #height constant
        bG=self.MHzScaleArr[np.argmax(self.imageAvgArr)]  #x axis offset
        cG=np.average(self.imageAvgArr[int(self.imageAvgArr.shape[0]/20.0):]) #Set the y axisoffsest guess based on the average of the last 10% of data
        dF0G=30 #good guess


        #index            0      1                2      3
        guess = np.array([aG    ,bG            , cG*1.0, dF0G ]) #array of guess



        PF, pcov = spo.curve_fit(self.spectralProfile, self.MHzScaleArr, self.imageAvgArr, p0=guess,sigma=error)
        self.F0=self.find_Voigt_Peak(self.MHzScaleArr,PF)
        dF0=PF[3] #FWHM of spectral profile

        def temp(freq):
            #freq: HHz
            return self.spectralProfile(freq,*PF)
        self.fitFunc=temp
        self.T=self.calculateTemp(dF0)


#---------SINGLE VOIGT FIT---------
#this voigt is normalized to a height of 1, then multiplied by the variable a
#there is no yaxis offset in the voigt, that is done later
    def voigt(self,f,a,f0,dF0,gamma):
        #units must be consistent!!
        #f, frequency value
        #a, height of voigt
        #f0, center frequency
        #dF0, FWHM of the gaussian
        #gamma, FWHM of the lorentzian
        sigma=dF0/(2*np.sqrt(2*np.log(2))) #convert gaussian FWHM to std
        gamma=gamma/2 #convert lorentzian FWHM to HWHM
        x=f-f0
        z=(x+gamma*1.0J)/(sigma*np.sqrt(2.0)) #complex argument
        V=np.real(spf.wofz(z))/(sigma*np.sqrt(2*np.pi)) #voigt profile


        #now normalize to 1 at peak, makes fitting easier
        z0=(gamma*1.0J)/(sigma*np.sqrt(2.0)) #z at peak
        V0=np.real(spf.wofz(z0))/(sigma*np.sqrt(2*np.pi)) #height of voigt at peak
        return a*V/V0 #makes the height equal to a


#------SEXTUPLE VOIGT FIT----------
    def spectralProfile(self,freq,a,b,c,dF0,gamma=gv.LiTau/1E6):
        #units must be consitant!
        #freq: frequency
        #a: height constant
        #b: center frequency
        #c: vertical offset
        #dF0: FWHM of doppler broadening
        #gamma: FWHM of lorentzian

        ratio=4*(1/.55)#ratio of intensity of f=2 to f=1. First number is power ratio in bands. Second fraction is ratio
                    #of hyperfine transitions
        a2=ratio*a
        a1=a
        val=c #start with the offset

        #F=2 transition
        val+=a2*self.voigt(freq, gv.S21, b+0*gv.F1Sep/1E6, dF0, gamma)
        val+=a2*self.voigt(freq, gv.S22, b+gv.F2Sep/1E6, dF0, gamma)
        val+=a2*self.voigt(freq, gv.S23, b+gv.F3Sep/1E6, dF0, gamma)


        #F=1 transition
        val+=a1*self.voigt(freq, gv.S10, b+gv.F0Sep/1E6, dF0, gamma)
        val+=a1*self.voigt(freq, gv.S11, b+gv.F1Sep/1E6, dF0, gamma)
        val+=a1*self.voigt(freq, gv.S12, b+gv.F2Sep/1E6, dF0, gamma)
        return val

    def calculateTemp(self,dF0,f0=gv.Li_D2_Freq,MHz=True):
        #calculate the temperatute in mk for a given doppler width at a given frequency F0
        if MHz==True:
            dF0=dF0*1E6
        return (dF0/f0)**2 *(gv.massLi7*gv.cLight**2)/(8*gv.kB*np.log(2))

    def cheatData(self,offset=5,showFit=True,returnFit=True):

        #need to include better interaction with the plot! doesn't work as expected when accessed from mainGUI


        smooth=spsi.savgol_filter(self.imageAvgArr,21,5)
        minHeight=(np.max(smooth)-np.min(smooth))/2.0 +np.min(smooth)
        peak=spsi.find_peaks(smooth,height=minHeight)[0][0]

        x=np.flip(self.MHzScaleArr[peak-offset:])
        y=np.flip(self.imageAvgArr[peak-offset:])
        #newPeak=184
        error=np.ones(x.shape[0])

        PF,stdDev=self.fit_Image_Data(x,y)
        a, b, c, dF0 = PF
        if showFit==True:
            voigtFit=[]
            for elem in x:
                voigtFit.append(self.spectralProfile(elem, a, b, c, dF0))

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
            voigtFit=np.asarray(self.spectralProfile(self.MHzScaleArr, a, b, c, dF0))
            return [a, b, c, dF0], stdDev, voigtFit
        return PF, stdDev

    def find_Voigt_Peak(self,MHzScaleArr,PF):
        x=np.linspace(self.MHzScaleArr[0],self.MHzScaleArr[-1],num=1000000)
        y=self.spectralProfile(x,*PF)
        F0=x[spsi.find_peaks(y)[0]]
        return F0[0]






