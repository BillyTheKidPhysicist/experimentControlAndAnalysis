import sys
import globalVariables as gv
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import fits
import scipy.special as sps
import scipy.optimize as spo
from MakeMHzScale import MHzScale
import scipy.ndimage as spni
from DataAnalysis import fitPixelData, spectralProfile
import scipy.interpolate as spi

def calculateTemp(dF0,f0=gv.Li_D2_Freq,MHz=True):
    #calculate the temperatute in mk for a given doppler FWHM width at a given frequency F0
    if MHz==True:
        dF0=dF0*1E6

    return (dF0/f0)**2 *(gv.massLi7*gv.cLight**2)/(8*gv.kB*np.log(2))


#----------------IMPORT AND RANGLE IMAGES INTO RIGHT SHAPE-----------
#Directory
path = "C:\Data\Runs\\2_9_21"
os.chdir(path)
fileName = 'run6Far'


#fileName='run15ShutterOpen'

#Opening fits file and creating array. Flip the y axis such that image is oriented correctly.
fitsFile = fits.open(fileName+'.fits')
imagesList = fitsFile[0].data
imagesArr=imagesList.astype(float)
imagesArr=np.flip(imagesArr,axis=1)




#--------------MAKE MHZ SCALE-------------------------

#get MHS array
infoFile=open(fileName[:-3]+'Info.txt')
startVolt=float(infoFile.readline().split(',')[1])
stopVolt=float(infoFile.readline().split(',')[1])
#print(startVolt,stopVolt)
dataFileName=fileName[:-3]+'DAQData.csv'
DAQData=np.loadtxt(dataFileName,delimiter=',')
MHZMaker=MHzScale(DAQData)
MHzScaleArr=MHZMaker.make_MHZ_Scale(False)

galvoVoltArr=DAQData[:, 0]
liRefVoltArr=DAQData[:, 1]
P=np.polyfit(galvoVoltArr, MHzScaleArr, 1)
#print(imagesArr.shape[0])
temp=np.linspace(startVolt, stopVolt, num=imagesArr.shape[0])
imageFreqMhzArr=P[1]+P[0]*temp


#----------IMAGE REFINEMENT--------------------
#choosing which images to add together, and what image region. This must be done carefully to capture enough information,
# but not add too much noise. also trimming hot and dead pixels




#There are two pixels right next to each other in the viewing region that gives us trouble. one is a hot pixel, the other
#is a cold pixel. Average them to the values of the surrounding pixels. ALso, clip the values of any other goofy pixels.
#after removing, smooth the image
for i in range(imagesArr.shape[0]):
    #this is only valid for a specific binning and field of view. Double check before using
    #pixHot=imagesArr[i][127,140]
    #pixCold=imagesArr[i][128,140]
    #pixHotAverage=(imagesArr[i][127,139]+imagesArr[i][127,141]+imagesArr[i][126,139]+imagesArr[i][126,141]+imagesArr[i][126,140])/5
    #pixColdAverage=(imagesArr[i][128,139]+imagesArr[i][128,141]+imagesArr[i][129,139]+imagesArr[i][129,141]+imagesArr[i][129,140])/5
    #imagesArr[i][127,140]=pixHotAverage
    #imagesArr[i][128,140]=pixColdAverage
    pass
    imagesArr[i]=spni.gaussian_filter(imagesArr[i], 0.0)


#this start and end value captures 90% of the beam in frequency space. It is used to make the single image that is the mean
#of a stack of iamges between imStart:imEnd
offset=-5
imStart=68+offset
imEnd=105+offset
imageMean=np.mean(imagesArr[imStart:imEnd],axis=0)







#crop the image to the focus to remove noise
xStart=85
xEnd=89
yStart=110
yEnd=130




#visualize the results


# #plot what the mean image looks like
temp=imageMean[yStart:yEnd,xStart:xEnd] #crop image to only the focus to remove noise

# plt.imshow(temp)
# plt.show()

# #plot the mean of the focus over the restricted images
temp=imagesArr[imStart:imEnd,yStart:yEnd,xStart:xEnd]
temp=np.mean(temp,axis=2)
temp=np.mean(temp,axis=1)
# plt.plot(temp)
# plt.grid()
# plt.show()


#--------------SPATIAL PROFILE OF LASER----------------------
#this can be done with either a voigt fit, or a Radial Basis Function (RBF) fit. The rbf is better for data
#that doesn't look very voigty

#expand the image region vertically
ya=yStart-30
yb=yEnd+30

magnification=2.8 #this can be found from using a ruler, or using optics rules.
pixelSize=magnification*5*.025 #in units of mm. binning X magnification X pixel size

#array of the profile
profArr=np.mean(imageMean[ya:yb,xStart:xEnd],axis=1) #laser profile from data


zArr=(np.arange(0,imagesArr[0].shape[0])+.5)*pixelSize #each position value corresponds to the center of the pixel




#two fit options, voigt or RBF

#------voigt fit
def voigtSpaceFit(x,x0,a,sigma,gamma,b):
    v0=sps.voigt_profile(0, sigma, gamma)
    v=sps.voigt_profile(x-x0, sigma, gamma)
    return a*v/v0+b
eps=1e-10 #small non zero number
guess=[55.0,50.0,5.0,5.0,1200]
bounds=[(-np.inf,eps,eps,eps,0.0),(np.inf,np.inf,np.inf,np.inf,np.inf)]
params, pcovLoopGammaFree=spo.curve_fit(voigtSpaceFit, zArr[ya:yb], profArr, p0=guess,bounds=bounds)
fwhm=.5346*(2*params[3])+np.sqrt(.2166*(2*params[3])**2+(params[2]*2.335)**2)
print('fwhm',fwhm,'Z0',params[0])
centerZVoigt=params[0]
plt.close('all')
plt.title('Vertical Spatial Profile Data With Voigt Fit')
plt.scatter(zArr[ya:yb],profArr,c='r',label='Data',s=50.0,marker='x')
plt.plot(zArr[ya:yb],profArr,c='r',linestyle=':')
plt.plot(zArr[ya:yb],voigtSpaceFit(zArr[ya:yb],*params),label='Fit')
plt.axvline(x=centerZVoigt,c='black',linestyle=':',label='Peak')
plt.grid()
plt.legend()
plt.xlabel('Position, mm')
plt.ylabel('Signal Strength, AU')
plt.savefig('Spatial_'+fileName)
# plt.show()


#----RBF fit-----
smoothFact=1.0#you might need to play with the smooth value. Change by factors of 10
rbfSpaceFit=spi.Rbf(zArr[ya:yb],profArr,smooth=smoothFact)
zArrFine=np.linspace(zArr[ya:yb].min(),zArr[ya:yb].max(),num=1000) #very fine z array
rbfFitArrFine=rbfSpaceFit(zArrFine) #very fine fit array
centerZRBF=zArrFine[np.argmax(rbfFitArrFine)]

#find fwhm from fit
maximum=np.max(rbfFitArrFine)

zLeftArg=np.argwhere(rbfFitArrFine>(maximum+np.min(rbfFitArrFine))/2)[0]
zRightArg=np.argwhere(rbfFitArrFine>(maximum+np.min(rbfFitArrFine))/2)[-1]

fwhm=zArrFine[zRightArg]-zArrFine[zLeftArg]
print('fwhm',fwhm,'Z0',centerZRBF)
plt.close('all')
plt.plot(zArr[ya:yb],profArr)
plt.plot(zArr[ya:yb],rbfSpaceFit(zArr[ya:yb]))
plt.axvline(x=centerZRBF,c='r',linestyle=':')
plt.axvline(x=zArrFine[zRightArg])
plt.axvline(x=zArrFine[zLeftArg])
plt.grid()
# plt.show()



#which center z to use
centerZ=centerZVoigt


#-------------FREQUENCY PROFILE AT EACH VERTICAL PIXEL------------------

def temperatureFitter(x, x0, a, sigma, gamma, b, singlePeak=True):
    #x: frequency, MHz
    #x0: center frequency, MHz
    #a: peak height, au
    #b: offset, au
    #sigma: standard deviation of gaussian contribution, MHz
    #gamma: HWHM of lorentzian, MHz
    #singlePeak: If true, assume optical pumping, if False use multipeak fit
    laserJitter=2.0/2.355
    sigma=np.sqrt(sigma**2+laserJitter**2)
    if singlePeak==True:
        v0=sps.voigt_profile(0,sigma,gamma)
        v=sps.voigt_profile(x-x0,sigma,gamma)
        return a*v/v0 +b
    else:
        val=spectralProfile(x,a,x0,b,dF0=sigma*2.355,gamma=2*gamma) #this function uses FWHM for linewidhts and constants are
        #ordered differently
        return val




valArr=np.mean(np.mean(imagesArr[:,yStart:yEnd,xStart:xEnd],axis=1),axis=1)
guess=[imageFreqMhzArr[np.argmax(valArr)], 100.0, 4.0, 5, 1200.0]
bounds=((-np.inf, 0, 0, 5.87/2.0, 0), (np.inf, np.inf, np.inf, 5.87*2, np.inf))
#print(valArr.shape)
params, pcov=spo.curve_fit(temperatureFitter, imageFreqMhzArr, valArr, p0=guess, bounds=bounds)
T=np.round(1e3*calculateTemp(2.355*params[2]),1)
#print( 'TEMP',T)

freqCenter=params[0]

plt.close('all')
plt.suptitle('Total Focus Signal vs Frequency')
plt.title('T= '+str(T)+' mk,   F0= '+str(np.round(freqCenter)) +' MHz')


xTest=np.linspace(imageFreqMhzArr[0], imageFreqMhzArr[-1], num=10000)

# plt.scatter(imageFreqMhzArr,valArr,c='r',label='Data',s=50.0,marker='x')
# plt.plot(imageFreqMhzArr,valArr,c='r',linestyle=':')
# plt.plot(xTest, temperatureFitter(xTest, *params),label='Fit')
# plt.axvline(x=freqCenter,c='black',linestyle=':',label='Peak')
# plt.xlabel('frequency, MHz')
# plt.grid()
# plt.legend()
# plt.xlabel('Frequency, MHz')
# plt.ylabel('Signal Strength, AU')
# plt.savefig('TotalTemperature_'+fileName)
# plt.show()


gammaList=[]
sigmaList=[]
TGammaFreeList=[]
TGammaLockedList=[]
temp_error = []
F0List=[]
yEnd+=3
yStart+=-3
for i in range(yStart,yEnd):
    #print(i)
    valArr=np.mean(imagesArr[:,i,xStart:xEnd],axis=1)
    #print(valArr.shape)
    def fitFreeGamma(x,x0,a,sigma,gamma,b):
        return temperatureFitter(x, x0, a, sigma, gamma, b)
    def fitLockedGamma(x,x0,a,sigma,b):
        gamma=5.87/2
        return temperatureFitter(x, x0, a, sigma, gamma, b)

    #-----gamma free--------
    guessGammaFree=[imageFreqMhzArr[np.argmax(valArr)], 100.0, 4.0, 5, 1200.0]
    boundsGammaFree=((-np.inf, 0, 0, 5.87/2.0, 0), (np.inf, np.inf, np.inf, np.inf, np.inf))
    paramsLoopGammaFree, pcovLoopGammaFree=spo.curve_fit(fitFreeGamma, imageFreqMhzArr, valArr, p0=guessGammaFree, bounds=boundsGammaFree)


    TGammaFree=1e3*calculateTemp(2.355*paramsLoopGammaFree[2])
    TGammaFreeList.append(TGammaFree)

    #-----gamma Locked--------
    guessGammaLocked=[imageFreqMhzArr[np.argmax(valArr)], 100.0, 4.0, 1200.0]
    boundsGammaLocked=((-np.inf, 0, 0, 0), (np.inf, np.inf, np.inf, np.inf))
    paramsLoopGammaLocked, pcovLoopGammaLocked=spo.curve_fit(fitLockedGamma, imageFreqMhzArr, valArr,
                                                             p0=guessGammaLocked, bounds=boundsGammaLocked)

    TGammaLocked=1e3*calculateTemp(2.355*paramsLoopGammaLocked[2])
    TGammaLockedList.append(TGammaLocked)

    #print(paramsLoopGammaFree[:])

    error=250*np.var(np.append(valArr[0:15],valArr[138:150]))
    print()
    print('Uncertainty',error**0.5)

    Temp = []
    chi_squared = []
    Coord = []
    SigmaUpperBound = []
    sigmaStart = 0.5
    sigmaEnd = 10
    sigmaIncrement = 0.05
    while sigmaStart <= sigmaEnd:

        def fitLockedSigma(x, x0, a, gamma, b):
            sigma = sigmaStart
            return temperatureFitter(x, x0, a, sigma, gamma, b)

        guessSigmaLocked=[imageFreqMhzArr[np.argmax(valArr)], 100.0, 8.0, 1200.0]
        boundsSigmaLocked=((-np.inf, 0, 0, 0), (np.inf, np.inf, np.inf, np.inf))
        paramsLoopSigmaLocked, pcovLoopSigmaLocked=spo.curve_fit(fitLockedSigma, imageFreqMhzArr, valArr,
                                                                 p0=guessSigmaLocked, bounds=boundsSigmaLocked)

        R2=[]

        for i in range(0, 150):
            value_fit=fitLockedSigma(imageFreqMhzArr[i], *paramsLoopSigmaLocked)
            value_data=valArr[i]
            difference=(value_data-value_fit)**2 / error
            #print(difference)
            R2.append(difference)

        chi_squared.append(sum(R2))
        #print('sigma',sigmaStart,'chi_squaured', chi_squared)
        Temp.append(1e3*calculateTemp(2.355*sigmaStart))
        Coord.append([sum(R2),1e3*calculateTemp(2.355*sigmaStart)])
        SigmaUpperBound.append(sigmaStart)
        sigmaStart += sigmaIncrement


    covTemp=1e3*calculateTemp(2.355*(paramsLoopGammaFree[2]+(pcovLoopSigmaLocked[2][2])**0.5))
    LowestTemp = Coord[np.argmin(Coord,axis=0)[0]][1]
    LowestChi = Coord[np.argmin(Coord,axis=0)[0]][0]


    def find_nearest(array, value):
        array=np.asarray(array)
        idx=(np.abs(array-value)).argmin()
        return array[idx]

    chi_plusone = find_nearest(chi_squared,LowestChi+1)
    chi_plusone_index = chi_squared.index(chi_plusone)
    sigmaTemp = Temp[chi_plusone_index]
    temp_error.append(np.abs(sigmaTemp-LowestTemp))
    sigmaFitUpperBound = SigmaUpperBound[chi_plusone_index]

    print('T From Chi', LowestTemp)
    print('T From Fit', TGammaFree)
    print('T Upper Value from Chi', sigmaTemp)
    print('T Upper Value from Cov', covTemp)
    #print(pcovLoopGammaFree[2][2])
    # plt.xlabel('Temperature, mK')
    # plt.ylabel('Chi squared')
    # plt.scatter(Temp,chi_squared)
    # plt.grid()
    # plt.show()

    #print('GammaLocked')
    #print(pcovLoopGammaLocked)
    #print(np.corrcoef(pcovLoopGammaLocked))

    # perr = np.sqrt(np.diag(pcovLoopGammaFree))
    # sigmaUpper = paramsLoopGammaFree[2] + 1.96*perr[2]
    # Tupper = 1e3*calculateTemp(2.355*sigmaUpper)
    # Tlower = 2*TGammaFree - Tupper
    # if Tlower < 0:
    #     Tlower = 0

    F0List.append(paramsLoopGammaFree[0])
    plt.close('all')
    plt.title('SINGLE ROW OF PIXELS IN FOCUS. Temp = %.2f mK' %sigmaTemp)
    plt.scatter(imageFreqMhzArr,valArr)
    xTest=np.linspace(imageFreqMhzArr[0],imageFreqMhzArr[-1],num=10000)
    plt.plot(xTest, temperatureFitter(xTest, paramsLoopGammaFree[0],paramsLoopGammaFree[1],sigmaFitUpperBound,paramsLoopGammaFree[3],paramsLoopGammaFree[4]))
    plt.grid()
    plt.xlabel('frequency, MHz')
    plt.show()


trim=10

F0Arr=np.asarray(F0List)
F0ArrTemp=F0Arr.copy()

zArrTemp=zArr[yStart:yEnd].copy()
F0ArrTemp=F0ArrTemp[trim:-trim] #trim out the 'bad' stuff
zArrTemp=zArrTemp[trim:-trim] #trim out the 'bad' stuff

vArrTemp=gv.cLight*F0ArrTemp*1e6/gv.Li_D2_Freq
zArrTemp=np.flip(zArrTemp) #because zarr starts at top of image, which is confusing

#since I am in a linear portion of the phase space plot I should be able to center things so absolute values
#don't matter
vArrTemp-=vArrTemp.mean() #center the velocity array at zero
zArrTemp-=np.mean(zArrTemp) #center the position array at zero
vArrTemp=1e3*vArrTemp #convert to mm/s


#fit the data with a line to find how much it must rotate by and thus the distance
#velocity = slope * distance, and thus 1/slope is how much time is required to focus!
# def phaseSpaceFunc(x,m,b):
#     return m*x+b
# #trim=15
# params,pcov=spo.curve_fit(phaseSpaceFunc,zArrTemp,vArrTemp)
# perr = np.sqrt(np.diag(pcov))
# v0=210.0 #atom speed
# dz=-np.round(1e3*v0/params[0]) #convert to mm
# dzError=np.round(dz*perr[0]/params[0]) #standard error propogation
# #print(dz,dzError)
#
# print(params)
# zValArr=np.asarray([25,55,70,100]) #mm
# dzValArr=np.asarray([176,152,158,157])

# plt.close('all')
# plt.title('Phase Space Plot, line params = '+str(np.round(params)))
# plt.plot(zArrTemp,phaseSpaceFunc(zArrTemp,*params))
# plt.scatter(zArrTemp,vArrTemp)
# plt.grid()
# plt.xlabel('Position, mm')
# plt.ylabel('Velocity, mm/s')
# plt.axvline(x=np.mean(zArrTemp),c='black',linestyle=':')
# plt.savefig('PhaseSpace'+str(fileName))
# plt.show()





dVArr=3e8*F0Arr*1e6/gv.Li_D2_Freq
#plt.plot(zArr[yStart:yEnd],TList)
plt.title('Temperature vs Vertical Position')
plt.title('Laser jitter=2.0MHz RMS')
#plt.scatter(zArr[yStart:yEnd], TGammaLockedList, label='gamma=5.87MHz')
plt.scatter(zArr[yStart:yEnd], TGammaFreeList, label='gamma not locked')
plt.errorbar(zArr[yStart:yEnd],TGammaFreeList,yerr=temp_error, linestyle="None",barsabove=True,capsize=3)
plt.xlabel('Vertical position in atomic beam, mm')
plt.ylabel('Temperature, mk')
plt.axvline(x=centerZ,c='r',linestyle=':',label='Atom beam center')
plt.grid()
plt.legend()
plt.savefig('TemperatureVertical_'+fileName)
plt.show()
