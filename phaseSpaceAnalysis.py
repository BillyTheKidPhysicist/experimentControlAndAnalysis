import scipy.signal as spsig
import sys
from DataAnalysis import fit_Spectral_Data
import globalVariables as gv
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import fits
import scipy.special as sps
import scipy.optimize as spo
from MakeMHzScale import MHzScale
import scipy.ndimage as spni
import scipy.interpolate as spi

def voigtImageFit(x, x0, a, sigma, gamma, b):
    v0=sps.voigt_profile(0, sigma, gamma)
    v=sps.voigt_profile(x-x0, sigma, gamma)
    return a*v/v0+b
def calculateTemp(dF0,f0=gv.Li_D2_Freq,MHz=True):
    #calculate the temperatute in mk for a given doppler FWHM width at a given frequency F0
    if MHz==True:
        dF0=dF0*1E6

    return (dF0/f0)**2 *(gv.massLi7*gv.cLight**2)/(8*gv.kB*np.log(2))

def laserSignal(x,P=1):
    #get the laser signal for our far field setup.
    #x: Position, mm
    #P: power in arbitrary units
    x0=33.2 #position in image along x axis in mm. Bottom left of image is defined as (0,0)
    w=30.7 # waist in mm
    return P*np.exp(-2*((x-x0)/w)**2)


#Directory
path = "C:\Data\Runs\\8_26_21"
os.chdir(path)
fileName = 'run42Far'



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





#Find background
imageBackGround=(np.mean(imagesArr[-5:],axis=0)+np.mean(imagesArr[:5],axis=0))/2 #take images from beginning and end
#to get average of background noise

imagesArr=imagesArr-imageBackGround #extract background.




#go through each image and apply a median filter to remove hot/dead pixels and cosmic rays
for i in range(imagesArr.shape[0]):
    imagesArr[i]=spni.median_filter(imagesArr[i], size=3)


#automatically find the peak frame
fwhmFact=1.5 #number of factors of the FWHM to include on either side of the frame

xStart=50
xEnd=-50
yStart=50
yEnd=-50
temp=np.mean(imagesArr[:,yStart:yEnd,yStart:yEnd],axis=2)
temp=np.mean(temp,axis=1)
temp=spsig.savgol_filter(temp,11,2)


guess=[np.argmax(temp),np.max(temp)-np.min(temp),1,1,np.min(temp)]
x=np.arange(temp.shape[0])
params,pcov=spo.curve_fit(voigtImageFit, x, temp, p0=guess)
fwhmImage=.5346*(2*params[3])+np.sqrt(.2166*(2*params[3])**2+(params[2]*2.335)**2)

#select the left and right image to capture most of the signal
centerImage=params[0]
imageLeft=int(centerImage-fwhmFact*fwhmImage)
imageRight=int(centerImage+fwhmFact*fwhmImage)
plt.axvline(x=imageLeft,c='black')
plt.axvline(x=imageRight,c='black')
plt.scatter(x,temp,c='r')
plt.plot(voigtImageFit(x, *params))
plt.grid()
plt.show()



#Spatial Analysis Box
yStart=0
yEnd=-1
delta=-50*0
xStart=110#109-delta
xEnd=130#111-delta
image=np.mean(imagesArr[imageLeft:imageRight,yStart:yEnd,xStart:xEnd],axis=0)
plt.imshow(image)
plt.show()


valArr=np.mean(np.mean(imagesArr[:,yStart:yEnd,xStart:xEnd],axis=1),axis=1)


lensHeating=False
laserJitter=1.1
spectralFit=fit_Spectral_Data(imageFreqMhzArr,valArr,lensHeating=lensHeating,vTMaxLens=.07*200*.9
                                         ,peakMode='single',laserJitter=laserJitter)


T=1e3*spectralFit.get_Temperature()
spectralFit.print_Results()


area = np.trapz(valArr,x=imageFreqMhzArr*10**6*2*np.pi) #4649

maxSignal=np.max(valArr)
print('area under curve in frequency space',area/1e6)
print('S(w0)',maxSignal/area)

freqCenter=spectralFit.fitResultsDict['center frequency']
print('Temperature of whole image, mk', np.round(T,1))

plt.close('all')
plt.suptitle('Total Focus Signal vs Frequency. Geometric broadening: '+str(lensHeating))
plt.title('T= '+str(np.round(T,1))+' mk,   F0= '+str(np.round(freqCenter)) +' MHz')
xTest=np.linspace(imageFreqMhzArr[0], imageFreqMhzArr[-1], num=10000)

#plt.scatter(imageFreqMhzArr,valArr,c='r',label='Data',s=50.0,marker='x')
plt.plot(imageFreqMhzArr,valArr,c='r')
plt.plot(imageFreqMhzArr, spectralFit.fit_Result_Function(imageFreqMhzArr),label='Fit')
plt.axvline(x=freqCenter,c='black',linestyle=':',label='Peak')
plt.xlabel('frequency, MHz')
plt.grid()
plt.legend()
plt.xlabel('Frequency, MHz')
plt.ylabel('Signal Strength, AU')
#plt.savefig('TotalTemperature_'+fileName)
plt.show()










Type='mean'
if Type=='sum':
    imageMean=np.sum(imagesArr[imageLeft:imageRight,yStart:yEnd,xStart:xEnd], axis=0)  #the stacked image
elif Type=='mean':
    imageMean=np.mean(imagesArr[imageLeft:imageRight,yStart:yEnd,xStart:xEnd], axis=0)  #the stacked image
else:
    raise Exception('INVALID')
print('Final stacked image is '+Type+' of images')




magnification=3  #this can be found from using a ruler, or using optics rules.
pixelSize=magnification*4*.024  #in units of mm. binning X magnification X pixel size

#array of the profile
profArr=np.mean(imageMean, axis=1)  #laser profile from data
zArr=(np.arange(0, profArr.shape[0]))*pixelSize  #each position value corresponds to the center of the pixel

def voigtSpaceFit(x, x0, a, sigma, gamma, b):
    v0=sps.voigt_profile(0, sigma, gamma)
    v=sps.voigt_profile(x-x0, sigma, gamma)
    return a*v/v0+b


eps=1e-10  #small non zero number
guess=[zArr[np.argmax(profArr)], profArr.max()-profArr.min(), 1.0, 1.0, profArr.min()]
print(guess)
bounds=[(-np.inf, eps, eps, eps, -np.inf), (np.inf, np.inf, np.inf, np.inf, np.inf)]
params=spo.curve_fit(voigtSpaceFit, zArr, profArr, p0=guess, bounds=bounds )[0]
x0,a,sigma,gamma,b=params
fwhm=.5346*(2*params[3])+np.sqrt(.2166*(2*params[3])**2+(params[2]*2.335)**2)

print('params',params)
print('fwhm', np.round(fwhm, 2), 'Z0', np.round(params[0],  1), 'sig peak', np.round(params[1]))
rMax=30 #bounds of integral in mm
rArr=np.linspace(0,rMax,num=1000)
# plt.plot(rArr,2*np.pi*rArr*voigtSpaceFit(rArr,0.0,1.0,sigma,gamma,0.0))
# plt.show()
# print("Integral",np.trapz(2*np.pi*rArr*voigtSpaceFit(rArr,0.0,1.0,sigma,gamma,0.0),x=rArr))
plt.title('Signal vs. Transerve Position')
plt.xlabel("Transverse position, mm")
plt.ylabel("Signal, au")
plt.scatter(zArr,profArr,c='r',label='Data')
plt.plot(zArr,voigtSpaceFit(zArr,*params),label='fit')
plt.legend()
plt.grid()
plt.show()











#
#
#
# #get the stacked image over the correct region
FWHM = []
FWHMSig = []
FWHMpos = []
width=8
i = 32
while i <= (256-32-width):

    xStart=0+i
    xEnd=xStart+width

    yStart=0
    yEnd=-1
    temp=imagesArr[0:-1,yStart:yEnd,xStart:xEnd]



    print('sum of images: ',np.round(np.sum(temp)/1e6,1))

    Type='mean'
    if Type=='sum':
        imageMean=np.sum(temp,axis=0) #the stacked image
    elif Type=='mean':
        imageMean=np.mean(temp,axis=0) #the stacked image
    else:
        raise Exception('INVALID')
    print('Final image is '+Type+' of images')

    # plt.imshow(imageMean)
    # plt.title('imageMean')
    # plt.show()

#
#
#
#
#
#
#     #--------------SPATIAL PROFILE OF LASER----------------------
#     #this can be done with either a voigt fit, or a Radial Basis Function (RBF) fit. The rbf is better for data
#     #that doesn't look very voigty
#
#
    magnification=3 #this can be found from using a ruler, or using optics rules.
    pixelSize=magnification*4*.024 #in units of mm. binning X magnification X pixel size

    #array of the profile
    profArr=np.mean(imageMean[yStart:yEnd],axis=1) #laser profile from data
    zArr=(np.arange(0,profArr.shape[0]))*pixelSize #each position value corresponds to the center of the pixel



    # def qGauss(z,z0,a,b,q,sigma):
    #
    #     eq=(1+(1-q)*(-((z-z0)/sigma)**2))**(1/(1-q))
    #     eq0=1#(1+(1-q)*(-sigma*(1e-10)**2))**(1/(1-q))
    #     return a*eq/eq0+b
    #
    # guess=[30,2000.0,5,1.00001,5]
    # params=spo.curve_fit(qGauss,zArr,profArr,p0=guess)[0]
    #
    # plt.title('Fit to data')
    # plt.xlabel("Transverse position, mm")
    # plt.ylabel("Signal, au")
    # plt.scatter(zArr,profArr,c='r',label='Data')
    # plt.plot(zArr,qGauss(zArr,*params),label='q-Gaussian fit')
    # plt.legend()
    # plt.grid()
    # plt.show()










    #------voigt fit
    def voigtSpaceFit(x,x0,a,sigma,gamma,b):
        v0=sps.voigt_profile(0, sigma, gamma)
        v=sps.voigt_profile(x-x0, sigma, gamma)
        return a*v/v0+b
    eps=1e-10 #small non zero number
    guess=[30.0,10000.0,1e-1,2.5,1240]
    bounds=[(-np.inf,eps,eps,eps,0.0),(np.inf,np.inf,np.inf,np.inf,np.inf)]
    params, pcovLoopGammaFree=spo.curve_fit(voigtSpaceFit, zArr, profArr, p0=guess,bounds=bounds,)
    fwhm=.5346*(2*params[3])+np.sqrt(.2166*(2*params[3])**2+(params[2]*2.335)**2)

    print('fwhm',np.round(fwhm,2),'Z0',np.round(params[0],1),'sig peak',np.round(params[1]))
    FWHM.append(fwhm)
    FWHMSig.append(params[1])
    FWHMpos.append((i-width/2)*0.29)
    i = i+width


fuckingarray = [FWHMpos,FWHM]

a_file=open("FWHMrun41.txt","w")
#np.savetxt(a_file,fuckingarray)

print(FWHMpos)
print(FWHM)
plt.scatter(FWHMpos,FWHM)
plt.title('FWHM of Voigt vs Position for '+str(fileName))
plt.xlabel("Distance from Face of Output of Magnet (mm)")
plt.ylabel("FWHM of Voigt (mm)")
#plt.savefig('FWHM_'+fileName)
plt.grid()
plt.show()

plt.scatter(FWHMpos,FWHMSig)
plt.title('Peak Signal of Voigt vs Position for '+str(fileName))
plt.xlabel("Distance from Left of Image (mm)")
plt.ylabel("Peak Signal of Voigt (arb.)")
plt.grid()
#plt.savefig('FWHMSignal_'+fileName)
plt.show()

# centerZVoigt=params[0]
# plt.close('all')
# plt.title('Vertical Spatial Profile Data With Voigt Fit')
# plt.scatter(zArr,profArr,c='r',label='Data')
# plt.plot(zArr,profArr,c='r',linestyle=':')
# plt.plot(zArr,voigtSpaceFit(zArr,*params),label='Fit')
# plt.axvline(x=centerZVoigt,c='black',linestyle=':',label='Peak')
# plt.grid()
# plt.legend()
# plt.xlabel('Position, mm')
# plt.ylabel('Signal Strength, AU')
# # plt.savefig('Spatial_'+fileName)
# plt.show()
#
#
# #which center z to use
# centerZ=centerZVoigt
#
#
# #-------------FREQUENCY PROFILE AND FIT OF ENTIRE IMAGE------------
# #the transverse velocity maximum for the lens is crucial to get right here
#
#
# valArr=np.mean(np.mean(imagesArr[:,yStart:yEnd,xStart:xEnd],axis=1),axis=1)
#
# dataAnalyzer=DataAnalyzer()
# params=dataAnalyzer.fit_Spectral_Profile(imageFreqMhzArr,valArr,lensHeating=True,vTMaxLens=8.5,peakMode='multi')
# T=1e3*dataAnalyzer.T #convert to mk
#
#
# freqCenter=params[0]
# print('Temperature of whole image, mk', np.round(T,1))
#
# plt.close('all')
# plt.suptitle('Total Focus Signal vs Frequency')
# plt.title('T= '+str(np.round(T,2))+' mk,   F0= '+str(np.round(freqCenter)) +' MHz')
# xTest=np.linspace(imageFreqMhzArr[0], imageFreqMhzArr[-1], num=10000)
#
# #plt.scatter(imageFreqMhzArr,valArr,c='r',label='Data',s=50.0,marker='x')
# plt.plot(imageFreqMhzArr,valArr,c='r')
# plt.plot(xTest, dataAnalyzer.spectral_Fit(xTest),label='Fit')
# plt.axvline(x=freqCenter,c='black',linestyle=':',label='Peak')
# plt.xlabel('frequency, MHz')
# plt.grid()
# plt.legend()
# plt.xlabel('Frequency, MHz')
# plt.ylabel('Signal Strength, AU')
# #plt.savefig('TotalTemperature_'+fileName)
# plt.show()



#FREQUENCY PROFILE AND FIT OF EACH COLUMN
'''

gammaList=[]
sigmaList=[]
TGammaFreeList=[]
TGammaLockedList=[]
F0List=[]
yStartLoop=50
yEndLoop=150
clumpPixels=2



for i in range(yStartLoop,yEndLoop):
    # print(i)
    valArr=np.mean(imagesArr[:,i:i+clumpPixels,xStart:xEnd],axis=1)
    valArr=np.mean(valArr,axis=1)
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





    F0List.append(paramsLoopGammaFree[0])

    # plt.close('all')
    # print('T',TGammaFree)
    # plt.title('SINGLE ROW OF PIXELS IN FOCUS')
    # plt.plot(imageFreqMhzArr,valArr)
    # xTest=np.linspace(imageFreqMhzArr[0],imageFreqMhzArr[-1],num=10000)
    # plt.plot(xTest, temperatureFitter(xTest, *paramsLoopGammaFree))
    # plt.grid()
    # plt.xlabel('frequency, MHz')
    # plt.show()


# trim=10
#
# F0Arr=np.asarray(F0List)
# F0ArrTemp=F0Arr.copy()
#
# zArrTemp=zArr[yStart:yEnd].copy()
# F0ArrTemp=F0ArrTemp[trim:-trim] #trim out the 'bad' stuff
# zArrTemp=zArrTemp[trim:-trim] #trim out the 'bad' stuff
#
# vArrTemp=gv.cLight*F0ArrTemp*1e6/gv.Li_D2_Freq
# zArrTemp=np.flip(zArrTemp) #because zarr starts at top of image, which is confusing
#
# #since I am in a linear portion of the phase space plot I should be able to center things so absolute values
# #don't matter
# vArrTemp-=vArrTemp.mean() #center the velocity array at zero
# zArrTemp-=np.mean(zArrTemp) #center the position array at zero
# vArrTemp=1e3*vArrTemp #convert to mm/s
#
#
# #fit the data with a line to find how much it must rotate by and thus the distance
# #velocity = slope * distance, and thus 1/slope is how much time is required to focus!
# def phaseSpaceFunc(x,m,b):
#     return m*x+b
# #trim=15
# params,pcov=spo.curve_fit(phaseSpaceFunc,zArrTemp,vArrTemp)
# perr = np.sqrt(np.diag(pcov))
# v0=205.0 #atom speed
# dz=-np.round(1e3*v0/params[0]) #convert to mm
# dzError=np.round(dz*perr[0]/params[0]) #standard error propogation
# #print(dz,dzError)
#
#
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





# dVArr=3e8*F0Arr*1e6/gv.Li_D2_Freq
#plt.plot(zArr[yStart:yEnd],TList)
plt.title('Temperature vs Vertical Position')
plt.title('Laser jitter=2.0MHz RMS')
plt.plot(zArr[yStartLoop:yEndLoop], TGammaLockedList, label='gamma=5.87MHz')
plt.plot(zArr[yStartLoop:yEndLoop], TGammaFreeList, label='gamma>=5.87MHz')
plt.xlabel('Vertical position in atomic beam, mm')
plt.ylabel('Temperature, mk')
plt.axvline(x=centerZ,c='r',linestyle=':',label='Atom beam center')
plt.grid()
plt.legend()
#plt.savefig('TemperatureVertical_'+fileName)
plt.show()



'''

