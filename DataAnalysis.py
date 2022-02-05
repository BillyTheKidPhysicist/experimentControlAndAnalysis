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
import warnings
import sigfig

#FUTURE PHYSICISTS, THIS USES OBJECT ORIENTED PROGRAMMING.
#Written by William Debenham(Huntington), billydebenham@gmail.com, 6197870499. Feel free to contact

class fitSolutionClass:
    def __init__(self,dataAnalyzerObject):
        self.fitResultsDict={} #dictionary to hold parameters from solution. This makes it much harder to get confused about
        #which value goes with which param.
        self._fit_params=dataAnalyzerObject.params
        self._unitsListForParamters=['MHz','au','au','MHz','MHz','MHz','','','m/s','mK']
        self._dataAnalyzerObject=dataAnalyzerObject
        self._fill_fitResultsDict()
        self.imageFreqArr=dataAnalyzerObject.imageFreqArr
        self.imageSignalArr=dataAnalyzerObject.imageSignalArr
    def get_Temperature(self):
        return self.fitResultsDict['temperature']
    def get_Params(self):
        return self._fit_params
        # pass
    def _fill_fitResultsDict(self):
        v0, a, b, sigma, gamma=self._dataAnalyzerObject.params

        self.fitResultsDict['center frequency']=v0
        self.fitResultsDict['signal height']=a
        self.fitResultsDict['vertical offset']=b
        self.fitResultsDict['sigma, STDev']=sigma
        self.fitResultsDict['gamma, FWHM']=gamma
        self.fitResultsDict['laser jitter']=self._dataAnalyzerObject.laserJitter
        self.fitResultsDict['peak mode']=self._dataAnalyzerObject.peakMode
        self.fitResultsDict['including lens heating']=False if self._dataAnalyzerObject.freqTMaxLens is None else True
        self.fitResultsDict['max transverse lens velocity']=self._dataAnalyzerObject.vTMaxLens
        self.fitResultsDict['temperature']=self._dataAnalyzerObject.T
    def print_Results(self):
        print(self)
    def fit_Function(self,freq,params):
        #this method allows you to pass your own params to the fit function. You can get the params from the fit,
        #then modify them as you desire. For example
        # params=fitSolution.get_Params()
        # params[1]=0.0
        # fit=fitSolution.(freqArr,*params)
        freqTMaxLens=self._dataAnalyzerObject.freqTMaxLens
        peakMode=self._dataAnalyzerObject.peakMode
        laserJitter=self._dataAnalyzerObject.laserJitter
        return self._dataAnalyzerObject._spectral_Profile(freq, *params, freqTMaxLens, peakMode, laserJitter)
    def fit_Result_Function(self,freq):
        return self._dataAnalyzerObject.spectral_Fit(freq)
    def __str__(self):
        #this is used when passing an instance of fitSolutionClass to print.
        precision=3
        assert len(self.fitResultsDict)==len(self._unitsListForParamters)
        stringToPrint='|----Spectral Fit Solution--------\n'
        stringToPrint+='|Solution Parameters (rounded to '+str(precision)+' sig figs): \n'
        i=0
        for name,value in self.fitResultsDict.items():
            print(name,value)
            unit=self._unitsListForParamters[i]
            if type(value)==bool or type(value)==str:
                if value==True and type(value) is not str:
                    stringValue='True'
                elif value==False and type(value) is not str:
                    stringValue='False'
                else:
                    stringValue=value
            elif value is None:
                stringValue='None'
            else:
                stringValue=str(sigfig.round(value,precision))
            stringToPrint+='|'+name+": "+stringValue+' '+unit+'\n'
            i+=1
        stringToPrint+='|---------------------------------'
        return stringToPrint

def fit_Spectral_Data(imageFreqArr, imageSignalArr, peakMode='multi', vTMaxLens=None
                             ,gammaMin=gv.LiTau/1E6,laserJitter=0.0,gammaLockedToNatural=False):

    dataAnalyzerObject=_DataAnalyzer()
    dataAnalyzerObject.fit_Spectral_Profile(imageFreqArr,imageSignalArr,peakMode,vTMaxLens,gammaMin,
                                                                                laserJitter,gammaLockedToNatural)
    fitSolutionObject=fitSolutionClass(dataAnalyzerObject)
    return fitSolutionObject

class _DataAnalyzer:
    def __init__(self):
        self.params=None #params for the current solution of fit_Spectral_Profile. contains [v0,a,b,sigma,gamma]
        #where v0: center frequency, a: height, b: offset from zero, sigma: standard deviation of temperature
        #disitribution, gamma: FWHM of natural linewidth plus other lorentzian broadening effects
        self.vTMaxLens=None
        self.gammaMin=None#provided
        self.T=None #temperature, kelvin
        self.F0=None #center frequency of the profile. This is not necesarily the peak of the fit because of the
        #underlying hyperfine structure
        self.fitPeak=None #to keep track of peak value of fit. This is necesary because more or less sampling can cause
        #the peak to be bigger or smaller slightly. FOr example, the profile may be fitted with poor resolution, then
        #when the user attempts to use those fit values to scale the data up, the scaled data can be taller
        self.freqTMaxLens=None #What the maximum transverse frequency shift of the atomic beam is from the geometric
        #'heating' from the lens, MHz
        self.peakMode=None #wether using single or multi peak fitting. Multi peak utilizes the 6 transitions
        self.laserJitter=None #laser jitter, MHz
        self.imageFreqArr=None #1D array of frequency values corresponding to each image
        self.imageSignalArr=None #1D array of means or sums of image regions of interest
        self.gammaLockedToNatural=None
    def fit_Spectral_Profile(self,imageFreqArr, imageSignalArr, peakMode,vTMaxLens,gammaMin,laserJitter,
                             gammaLockedToNatural):
        #Fit spectral data. Can be (single peak or multipeak) and include the velocity distribution of lens.
        #imageFreqArr: frequency, or anything else, scale
        #imageSignalArr: Flourescence signal
        #peakMode: Use either multipeak (6 peaks) or single peak.
        #vTMaxLens: transverse velocity maximum at lens output. This is used to calculate the geoemtric 'heating'.
        #If None don't apply anything'
        #gammeMin: value of gamma. Set to None to enable free fitting of gamma. Set to a specific value to lock to that
        #value, or allow to go above depending on gammaFloor
        #laserJitter: Jitter of laser, standard deviation, MHz.
        #gammaLockedToNatural: set gamma to natural linewidth


        if vTMaxLens is None:
            self.freqTMaxLens=None
        else:
            self.freqTMaxLens=(vTMaxLens/gv.cLight)*gv.Li_D2_Freq/1e6 #transverse frequency maximum lens, MhZ
        gammaMax=30.0
        assert gammaMin<gammaMax
        if gammaLockedToNatural==True:
            assert gammaMin==gv.LiTau/1e6
        self.gammaLockedToNatural=gammaLockedToNatural
        self.peakMode=peakMode
        self.laserJitter=laserJitter
        self.vTMaxLens=vTMaxLens
        self.gammaMin=gammaMin
        self.imageFreqArr=imageFreqArr
        self.imageSignalArr=imageSignalArr
        self.catch_Errors_And_Give_Warnings()

        def minimize(params):
            #wrapper function to be minimized. One cannot do least squares for the curve fit because if the lensHeating
            #feature is enable, there will be a convolution, which has to be evaluated all at once for every data point.
            #using the scipy curve_fit this would be far to slow. It's better to compare entire fits at once, rather
            #than compare a fit point by point which requires recomputing the convolution over an dover again.
            v0, a, b, sigma, gamma=params
            fit=self._spectral_Profile(self.imageFreqArr,v0,a,b,sigma,gamma,self.freqTMaxLens,self.peakMode,self.laserJitter)
            return self._cost(imageSignalArr,fit)
        deltaFMax=100
        F0Guess=self.imageFreqArr[np.argmax(self.imageSignalArr)]
        aGuess=np.max(self.imageSignalArr)-np.min(self.imageSignalArr)
        bGuess=np.min(self.imageSignalArr)
        if self.gammaLockedToNatural==True:
            gammaMax=gammaMin+1e-6
        bounds=[(F0Guess-deltaFMax,F0Guess+deltaFMax),(0.0,aGuess*2),(bGuess-aGuess,bGuess+aGuess),(0,30),(gammaMin,gammaMax)]
        np.random.seed(42) #seed the generator to make results repeatable
        sol=spo.differential_evolution(minimize,bounds,polish=False) #a type of genetic algorithm, very robust.
        np.random.seed(int(time.time())) #resead to make pseudorandom
        self.params=sol.x
        self.T=self._calculate_Temp(self.params[3])
        self.F0=self.params[0]

        print(bounds)
        print(self.params)
        xTest=np.linspace(imageFreqArr[0],imageFreqArr[-1],num=1000)
        plt.grid()
        plt.plot(imageFreqArr,imageSignalArr,marker='x')
        plt.plot(xTest,self.spectral_Fit(xTest))
        plt.show()
    def catch_Errors_And_Give_Warnings(self):
        if self.peakMode!='multi' and self.peakMode!='single':
            raise Exception('Invalid peak mode. Choose \'multi\' or \'single\'')
        if not 100.0>self.gammaMin>=0.0:
            raise Exception('Gamma must be in MHz and positive and it appears not to be')
        if self.gammaMin<gv.LiTau/1e6:
            warnings.warn("Gamma is set to a value below the natural linewidth")
        if not 100.0>self.laserJitter>=0.0:
            raise Exception('Laser jitter must be in MHz and positive and it appears not to be')
        maxReasonableVTransLen=20.0
        if self.vTMaxLens is not None and (not maxReasonableVTransLen>self.vTMaxLens>=0.0):
            raise Exception('Transverse particle velocity leaving lens is set to a value that does not make sense. \n'
                            'must be m/s, positive and less than a reasonably sized value')
        maxExpectedFreqInMHz=5000
        if np.max(np.abs(self.imageFreqArr))>1e6*maxExpectedFreqInMHz:
            raise Exception('It appears that the provided frequency data is not in MHz')
        if isinstance(self.imageFreqArr, np.ndarray)==False or isinstance(self.imageSignalArr, np.ndarray)==False:
            raise Exception('Frequency data and signal data must be numpy arrays')
    def spectral_Fit(self,v):
        #the function that returns the fit to the data provided in self.fit_Spectral_Profile.
        #v: frequency, MHz
        #return: The "signal", in arbitrary units
        return self._spectral_Profile(v,*self.params,self.freqTMaxLens,self.peakMode,self.laserJitter)
    def _spectral_Profile(self,v,v0,a,b,sigma,gamma,freqTMaxLens,peakMode,laserJitter):
        #todo: this can be made more clean and logical with lens stuff
        #private method that returns the spectral profile. It may be convoluted with the geometric "heating" and or
        #contain multiple peaks

        profile=np.zeros(v.shape)+b
        sigma=np.sqrt(sigma**2+laserJitter**2) #gaussian standard deviation. laserJitter is assumed to be gaussian here

        if peakMode=='multi':
            profile+=self._multi_Voigt(v, v0, a, sigma, gamma = gamma)
        else:
            profile+=self._voigt(v, v0, a, sigma, gamma)
        if freqTMaxLens is not None:
            # plt.plot(v,profile)
            profileMax=profile.max() #to aid in normalizing and rescaling the profile after convolution
            profileMin=profile.min() #to aid in normalizing and rescaling the profile after convolution
            profile=(profile-profileMin)/(profileMax-profileMin)
            vLens=np.linspace(-(v.max()-v.min())/2,(v.max()-v.min())/2,num=v.shape[0]) #need to have the lens profile
            #centered for the convolution to preserve the center of the spectral profile that results
            lensVelProfile=np.vectorize(self.lens_Velocity_Spread)(vLens,freqTMaxLens) #not very efficient
            # plt.plot(v,lensVelProfile)
            profile=np.convolve(lensVelProfile,profile,mode='same') #convolution has distributive property so don't need

            profile=(profile-profile.min())/(profile.max()-profile.min())

            # plt.plot(v,profile)
            # plt.show()
            #to perform on each peak individually
            profile=profile*(profileMax-profileMin)+profileMin #rescale the peak
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
        #x0: maximum transverse velocity or frequency. Same units as x
        if np.abs(x)>x0:  #to prevent imaginary
            return 0
        else:
            return a*np.sqrt(1-(x/x0)**2)

        #---------SINGLE VOIGT FIT---------
        #this voigt is normalized to a height of 1, then multiplied by the variable a
        #there is no yaxis offset in the voigt, that is done later
    def _voigt(self,f, f0, a, sigma, gamma):
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

    def _multi_Voigt(self, freq, freq0, a, sigma, gamma=gv.LiTau/1E6):
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
        val=0.0


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
        if MHz==True:
            f0=f0/1e6

        sigmaVelocity=(sigma/f0)*gv.cLight
        toMkelvin=1e3
        return toMkelvin*sigmaVelocity**2*gv.massLi7/gv.kB
    def calculate_Temp(self,sigma,MHz=True):
        return self._calculate_Temp(sigma,MHz=MHz)

    def self_Test(self):

        #_spectral_Profile(self,v,v0,a,b,sigma,gamma,vTMaxLens,peakMode,laserJitter)
        numTestDataPoints=100
        testDataFreq=np.linspace(-100,100,num=numTestDataPoints)
        freq0=10
        sigma=5.0
        gamma=10.0
        b=1000.0
        a=500.0
        # testDataSignal=500.0*np.exp(-testDataFreq**2/100)+1000.0
        testDataSignal=spf.voigt_profile(testDataFreq-freq0,sigma,gamma)/spf.voigt_profile(0.0,sigma,gamma)
        testDataSignal=a*testDataSignal+b
        np.random.seed(42) #set the seed to repeatable noise

        # testDataSignal+=30.0*np.random.random(numTestDataPoints) #add gaussion noise
        np.random.seed(int(time.time())) #reseed the random generator

        tester1=fit_Spectral_Data(testDataFreq,testDataSignal,'multi',False,10.0,gv.LiTau/1e6,0.0)._dataAnalyzerObject
        tester2=fit_Spectral_Data(testDataFreq,testDataSignal,'single',False,10.0,gv.LiTau/1e6,0.0)._dataAnalyzerObject
        tester3=fit_Spectral_Data(testDataFreq,testDataSignal,'multi',True,10.0,gv.LiTau/1e6,0.0)._dataAnalyzerObject
        tester4=fit_Spectral_Data(testDataFreq,testDataSignal,'single',True,10.0,gv.LiTau/1e6,0.0)._dataAnalyzerObject
        tester5=fit_Spectral_Data(testDataFreq,testDataSignal,'single',True,10.0,gv.LiTau/1e6,1.0)._dataAnalyzerObject
        testerList=[tester1,tester2,tester3,tester4,tester5]
        # now test that saved results match with current results
        testResults=np.loadtxt('DataAnalyzer_TestData')
        for i in range(len(testerList)):
            assert np.all(testResults[i]==np.asarray([*testerList[i].params, testerList[i].T]))
            print('test '+str(i)+' passed')
        print('self test passed')
    def save_Test_Data(self,testerList):
        gv.warning_Sound()
        # input("you are trying to save new test data, and overwrite the old. Press enter to proceed")
        data=[]
        for tester in testerList:
            data.append([*tester.params,tester.T])
        np.savetxt('DataAnalyzer_TestData',np.asarray(data))

def perform_Test():
    tester=_DataAnalyzer()
    tester.self_Test()

# perform_Test()
