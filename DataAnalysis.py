import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as spo
import scipy.special as spf
import scipy.interpolate as spi
import scipy.signal as spsi
import globalVariables as gv
import time
import numba

#FUTURE PHYSICISTS, THIS USES OBJECT ORIENTED PROGRAMMING.
#Written by William Debenham(Huntington), billydebenham@gmail.com, 6197870499. Feel free to contact


#TODO: CHANGE FITFUNC TO IT'S OWN COMMENTED FUNCTION


class DataAnalyzer:
    def __init__(self):
        self.params=None #params for the current solution of fit_Spectral_Profile. contains [v0,a,b,sigma,gamma]
        #where v0: center frequency, a: height, b: offset from zero, sigma: standard deviation of temperature
        #disitribution, gamma: FWHM of natural linewidth plus other lorentzian broadening effects



        self.T=None #temperature, kelvin
        self.F0=None #center frequency of the profile
        self.fitPeak=None #to keep track of peak value of fit. This is necesary because more or less sampling can cause
        #the peak to be bigger or smaller slightly. FOr example, the profile may be fitted with poor resolution, then
        #when the user attempts to use those fit values to scale the data up, the scaled data can be taller
        self.freqTMaxLens=None #What the maximum transverse frequency shift of the atomic beam is from the geometric
        #'heating' from the lens, MHz
        self.peakMode=None #wether using single or multi peak fitting. Multi peak utilizes the 6 transitions
        self.laserJitter=None #laser jitter, MHz
        self.imageFreqArr=None #1D array of frequency values corresponding to each image
        self.imageSignalArr=None #1D array of means or sums of image regions of interest
    def fit_Spectral_Profile(self,imageFreqArr, imageSignalArr, peakMode='multi', lensHeating=False,vTMaxLens=10.0
                             ,gamma=gv.LiTau/1E6,laserJitter=0.0):
        #Fit spectral data. Can be (single peak or multipeak) and include the velocity distribution of lens.
        #imageFreqArr: frequency, or anything else, scale
        #imageSignalArr: Flourescence signal
        #peakMode: Use either multipeak (6 peaks) or single peak.
        #lensHeating: Include the effects of the convolution with the lens output transverse velocity distribution. This
        #is a form of fake 'heating', so should be deconvolved from the signal.
        #vTMaxLens: transverse velocity maximum at lens output. This is used to calculate the geoemtric 'heating'
        #gamme: value of gamma. Set to None to enable free fitting of gamma. Set to a specific value to lock to that
        #value, or allow to go above depending on gammaFloor
        #laserJitter: Jitter of laser, standard deviation, MHz.

        if peakMode!= 'multi' and peakMode!='single':
            raise Exception('Invalid peak mode. Choose multipeak or single peak')
        if lensHeating==False:
            freqTMaxLens=None
        else:
            freqTMaxLens=(vTMaxLens/gv.cLight)*gv.Li_D2_Freq/1e6 #transverse frequency maximum lens, MhZ
        self.fitPeak=None #reset the fit peak value
        self.freqTMaxLens=freqTMaxLens
        self.peakMode=peakMode
        self.laserJitter=laserJitter

        self.imageFreqArr=imageFreqArr
        self.imageSignalArr=imageSignalArr


        def minimize(params):
            #wrapper function to be minimized. One cannot do least squares for the curve fit because if the lensHeating
            #feature is enable, there will be a convolution, which has to be evaluated all at once for every data point.
            #using the scipy curve_fit this would be far to slow. It's better to compare entire fits at once, rather
            #than compare a fit point by point which requires recomputing the convolution over an dover again.
            v0, a, b, sigma, gamma=params
            fit=self._spectral_Profile(self.imageFreqArr,v0,a,b,sigma,gamma,freqTMaxLens,peakMode,laserJitter)
            return self._cost(imageSignalArr,fit)
        deltaF=100
        F0Guess=self.imageFreqArr[np.argmax(self.imageSignalArr)]
        aGuess=np.max(self.imageSignalArr)-np.min(self.imageSignalArr)
        bGuess=np.min(self.imageSignalArr)

        bounds=[(F0Guess-deltaF,F0Guess+deltaF),(0.0,aGuess*2),(-2*np.abs(bGuess),2*np.abs(bGuess)),(0,30),(gamma,30)]

        np.random.seed(42) #seed the generator to make results repeatable
        sol=spo.differential_evolution(minimize,bounds,polish=False) #a type of genetic algorithm, very robust.

        np.random.seed(int(time.time())) #resead to make pseudorandom
        self.params=sol.x
        self.T=self._calculate_Temp(self.params[3])
        self.F0=self.params[0]
        self.fitPeak=self.spectral_Fit(imageFreqArr).max()

        # print(self.params)
        # print(bounds)
        # plt.grid()
        # plt.plot(imageFreqArr,imageSignalArr)
        # plt.plot(imageFreqArr,self.spectral_Fit(imageFreqArr))
        # plt.show()
        return self.params

    def spectral_Fit(self,v):
        #the function that returns the fit to the data provided in self.fit_Spectral_Profile.
        #v: frequency, MHz
        #return: The "signal", in arbitrary units

        if self.params is None:
            raise Exception("You have not fit any data yet. You must call fit_Spectral_Profile first")
        return self._spectral_Profile(v,*self.params,self.freqTMaxLens,self.peakMode,self.laserJitter)
    def _spectral_Profile(self,v,v0,a,b,sigma,gamma,vTMaxLens,peakMode,laserJitter):
        #private method that returns the spectral profile. It may be convoluted with the geometric "heating" and or
        #contain multiple peaks

        profile=np.zeros(v.shape)
        sigma=np.sqrt(sigma**2+laserJitter**2) #gaussian standard deviation. laserJitter is assumed to be gaussian here

        if peakMode=='multi':
            profile+=self._multi_Voigt(v, v0, a, b, sigma, gamma)
        else:
            profile+=self._voigt(v, v0, a, sigma, gamma)
        if vTMaxLens is not None:
            peak0=profile.max() #to aid in normalizing and rescaling the profile after convolution
            profile=profile/peak0
            vLens=np.linspace(-(v.max()-v.min())/2,(v.max()-v.min())/2,num=v.shape[0]) #need to have the lens profile
            #centered for the convolution to preserve the center of the spectral profile that results
            lensVelProfile=np.vectorize(self.lens_Velocity_Spread)(vLens,vTMaxLens) #not very efficient
            profile=np.convolve(lensVelProfile,profile,mode='same') #convolution has distributive property so don't need
            #to perform on each peak individually
            if self.fitPeak is None:
                profile=profile*peak0/profile.max() #rescale the peak
            else:
                profile=self.fitPeak*profile/profile.max()
        return profile
    @staticmethod
    @numba.njit(numba.float64(numba.float64[:],numba.float64[:]))#special compiler to make code run faster
    def _cost(imageSignalArr,fit):
        #cost function for the differential evolution analysis
        return np.sqrt(np.sum((imageSignalArr-fit)**2))
    @staticmethod
    @numba.njit(numba.float64(numba.float64,numba.float64,numba.float64,numba.float64))#special compiler to make code run faster
    def _gauss(T,v,m=gv.massLi7,v0=0.0):
        t1=np.sqrt(m/(2*np.pi*gv.kB*T))
        t2=np.exp(-m*np.power(v-v0,2)/(2*gv.kB*T))
        return t1*t2
    @staticmethod
    @numba.njit() #special compiler to make code run faster
    def lens_Velocity_Spread(x, x0, a=1):
        #1D transvers velocity distribution for the output of the lens. This is because the lens a circular input
        #and so creates a unique semicircular distirbution. Can be derived considering the density and y velocity as a
        # function of y in the circular aperture easily.
        #x: either velocity or frequency value. Must have same units as x0
        #x0: maximum transverse velocity of frequency. Same units as x
        if np.abs(x)>x0:  #to prevent imaginary
            return 0
        else:
            return a*np.sqrt(1-(x/x0)**2)

        #---------SINGLE VOIGT FIT---------
        #this voigt is normalized to a height of 1, then multiplied by the variable a
        #there is no yaxis offset in the voigt, that is done later
    @staticmethod
    def _voigt(f, f0, a, sigma, gamma):
        #units must be consistent!!
        #f, frequency value
        #f0, center frequency
        #a, height of voigt
        #sigma,standard deviation of the gaussian
        #gamma, FWHM of the lorentzian
        gamma=gamma/2  #convert lorentzian FWHM to HWHM
        x=f-f0
        z=(x+gamma*1.0J)/(sigma*np.sqrt(2.0))  #complex argument
        V=np.real(spf.wofz(z))/(sigma*np.sqrt(2*np.pi))  #voigt profile

        #now normalize to 1 at peak, makes fitting easier
        z0=(gamma*1.0J)/(sigma*np.sqrt(2.0))  #z at peak
        V0=np.real(spf.wofz(z0))/(sigma*np.sqrt(2*np.pi))  #height of voigt at peak
        return a*V/V0  #makes the height equal to a

        #------SEXTUPLE VOIGT FIT----------

    def _multi_Voigt(self, freq, freq0, a, b, sigma, gamma=gv.LiTau/1E6):
        #units must be consitant!
        #freq: frequency
        #a: height constant
        #b: center frequency
        #c: vertical offset
        #sigma: standard deviation of gaussian profile
        #gamma: FWHM of lorentzian
        aRatio=4*gv.F2F1Ratio  #ratio of intensity of f=2 to f=1. First number is power ratio in sideband. Second
        # fraction is ratio of hyperfine transitions (see globalVariable.py for more info in the comments for
        #F2F1Ratio).
        a2=(aRatio/(aRatio+1))*a #I do some funny business here to try and get the parameter "a" to more closely match
        #the total height. a2/a1 still equals the parameter "aRatio", but now they also add up the parameter "a"
        a1=a*1/(aRatio+1) #same funny business
        val=b  #start with the offset


        #F=2 transition
        val+=a2*self._voigt(freq, freq0+gv.F1Sep/1E6, gv.S21, sigma, gamma)
        val+=a2*self._voigt(freq, freq0+gv.F2Sep/1E6, gv.S22, sigma, gamma)
        val+=a2*self._voigt(freq, freq0+gv.F3Sep/1E6, gv.S23, sigma, gamma)
        #F=1 transition
        val+=a1*self._voigt(freq, freq0+gv.F0Sep/1E6, gv.S10, sigma, gamma)
        val+=a1*self._voigt(freq, freq0+gv.F1Sep/1E6, gv.S11, sigma, gamma)
        val+=a1*self._voigt(freq, freq0+gv.F2Sep/1E6, gv.S12, sigma, gamma)
        return val

    def _calculate_Temp(self,sigma,f0=gv.Li_D2_Freq, MHz=True):
        dF0=2*np.sqrt(2*np.log(2))*sigma
        if MHz==True:
            dF0=dF0*1E6
        return (dF0/f0)**2*(gv.massLi7*gv.cLight**2)/(8*gv.kB*np.log(2))
    def calculate_Temp(self,sigma,MHz=True):
        return self._calculate_Temp(sigma,MHz=MHz)


