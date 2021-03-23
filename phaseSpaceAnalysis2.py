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
from DataAnalysis import fitPixelData

def calculateTemp(dF0,f0=gv.Li_D2_Freq,MHz=True):
    #calculate the temperatute in mk for a given doppler FWHM width at a given frequency F0
    if MHz==True:
        dF0=dF0*1E6

    return (dF0/f0)**2 *(gv.massLi7*gv.cLight**2)/(8*gv.kB*np.log(2))


#----------------IMPORT AND RANGLE IMAGES INTO RIGHT SHAPE-----------



fileName='run15ShutterOpen'

#Opening fits file and creating array. Flip the y axis such that image is oriented correctly.
fitsFile = fits.open(fileName+'.fits')
imagesList = fitsFile[0].data
imagesArr=imagesList.astype(float)
imagesArr=np.flip(imagesArr,axis=1)





#
# def fit(x,x0,a,sigma,gamma,b):
#     laserJitter=3
#     gamma+=5.6/2.0
#     sigma=np.sqrt(sigma**2+laserJitter**2)
#     z0=sps.voigt_profile(0,sigma,gamma)
#     z=sps.voigt_profile(x-x0,sigma,gamma)
#     return a*z/z0 +b
#
#
#
# valArr=np.mean(np.mean(imagesArr[:,yStart:yEnd,xStart:xEnd],axis=1),axis=1)
#
# print(1e3*calculateTemp(fitPixelData(imageFreqMhzArr,valArr)[-1]))
#
#
#
#
#
#
# guess=[imageFreqMhzArr[np.argmax(valArr)], 100.0, 4.0, 5.0, 1140.0]
# bounds=((-np.inf, 0, 0, 0, 0), (np.inf, np.inf, np.inf, np.inf, np.inf))
#
# #print(valArr.shape)
# params, pcov=spo.curve_fit(fit, imageFreqMhzArr, valArr, p0=guess, bounds=bounds)
# print('F0',params[0],'TEMP',1e3*calculateTemp(2.355*params[2]))
# #params=guess
#
# plt.title('FULL IMAGE OF FOCUS')
# plt.plot(imageFreqMhzArr,valArr)
#
# xTest=np.linspace(imageFreqMhzArr[0], imageFreqMhzArr[-1], num=10000)
# plt.plot(xTest, fit(xTest, *params))
# plt.xlabel('frequency, MHz')
# plt.show()
#
# sys.exit()

#--------------MAKE MHZ SCALE-------------------------
imageFreqMhzArr=np.linspace(0,1.0,num=imagesArr.shape[0])*.3*565.0
#get MHS array


#----------IMAGE REFINEMENT--------------------
#choosing which images to add together, and what image region. This must be done carefully to capture enough information,
# but not add too much noise. also trimming hot and dead pixels






upperTrim=1.5e3
imagesArr[imagesArr>upperTrim]=upperTrim
lowerTrim=1.0e3
imagesArr[imagesArr<lowerTrim]=lowerTrim

imStart=0
imEnd=-1
xStart=270
xEnd=277
yStart=260
yEnd=300






#There are two pixels right next to each other in the viewing region that gives us trouble. one is a hot pixel, the other
#is a cold pixel. Average them to the values of the surrounding pixels. ALso, clip the values of any other goofy pixels.
#after removing, smooth the image
for i in range(imagesArr.shape[0]):
    imagesArr[i]=spni.gaussian_filter(imagesArr[i], 1)


imageMean=np.mean(imagesArr[imStart:imEnd],axis=0)#[yStart:yEnd,xStart:xEnd]
# plt.imshow(imageMean[yStart:yEnd,xStart:xEnd])
# plt.show()

#this start and end value captures 90% of the beam in frequency space. It is used to make the single image that is the mean
#of a stack of iamges between imStart:imEnd

imStart=15
imEnd=45




#visualize the results


# #
# # #plot the mean of the focus over the restricted images
temp=imagesArr[imStart:imEnd,yStart:yEnd,xStart:xEnd]
temp=np.mean(temp,axis=2)
temp=np.mean(temp,axis=1)
# plt.plot(temp)
# plt.grid()
# plt.show()


#--------------SPATIAL PROFILE OF LASER----------------------

ya=yStart-50
yb=yEnd+50
pixelSize=2*4*.025 #in units of mm. binning X magnification X pixel size

profData=np.mean(imageMean[ya:yb,xStart:xEnd],axis=1) #laser profile from data
zArr=(np.arange(0,imagesArr[0].shape[0])+.5)*pixelSize #each position value corresponds to the center of the pixel
#print(zArr)
# plt.plot(zArr[ya:yb],profData)
# plt.show()

def fit(x,x0,a,b,c):
    #simple gaussian fit
    return a*np.exp(-(x-x0)**2/(2.0*b**2))+c
guess=[55.0,50.0,5,1200]
params,pcov=spo.curve_fit(fit,zArr[ya:yb],profData,p0=guess)
# print('fwhm',params[2]*2.335,'Z0',params[0])
# print('a',params[1])
# print(np.sqrt(np.diag(pcov))[1])
centerZ=params[0]
# plt.plot(zArr[ya:yb],profData)
# plt.plot(zArr[ya:yb],fit(zArr[ya:yb],*params))
# plt.grid()
# plt.show()


#-------------FREQUENCY PROFILE AT EACH VERTICAL PIXEL------------------

def fit(x,x0,a,sigma,gamma,b):
    laserJitter=2.0
    sigma=np.sqrt(sigma**2+laserJitter**2)
    v0=sps.voigt_profile(0,sigma,gamma)
    v=sps.voigt_profile(x-x0,sigma,gamma)
    return a*v/v0 +b


valArr=np.mean(np.mean(imagesArr[:,yStart:yEnd,xStart:xEnd],axis=1),axis=1)
guess=[imageFreqMhzArr[np.argmax(valArr)], 100.0, 4.0, 5.0, 1200.0]
bounds=((-np.inf, 0, 0, 0.0, 0), (np.inf, np.inf, np.inf, np.inf, np.inf))
#print(valArr.shape)
params, pcov=spo.curve_fit(fit, imageFreqMhzArr, valArr, p0=guess, bounds=bounds)
#print('F0',params[0],'TEMP',1e3*calculateTemp(2.355*params[2]))
#
# plt.title('FULL IMAGE OF FOCUS')
# plt.plot(imageFreqMhzArr,valArr)
#
# xTest=np.linspace(imageFreqMhzArr[0], imageFreqMhzArr[-1], num=10000)
# plt.plot(xTest, fit(xTest, *params))
# plt.xlabel('frequency, MHz')
# plt.show()


gammaList=[]
sigmaList=[]
TFreeList=[]
TFreeErrorList=[]
TLockedList=[]
TLockedErrorList=[]
F0List=[]
for i in range(yStart,yEnd):
    #print(i)
    valArr=np.mean(imagesArr[:,i,xStart:xEnd],axis=1)
    #print(valArr.shape)
    def fitFreeGamma(x,x0,a,sigma,gamma,b):
        return fit(x,x0,a,sigma,gamma,b)
    def fitLockedGamma(x,x0,a,sigma,gamma,b):
        gamma=5.87/2.0
        return fit(x,x0,a,sigma,gamma,b)
    params,pcov=spo.curve_fit(fitFreeGamma,imageFreqMhzArr,valArr,p0=guess,bounds=bounds)
    perr=np.sqrt(np.diag(pcov))
    T=1e3*calculateTemp(2.355*params[2])
    print(params[-3])
    TUpper=1e3*calculateTemp(2.355*(params[2]+perr[2]))
    TLower=1e3*calculateTemp(2.355*(params[2]-perr[2]))
    TFreeList.append(T)
    TFreeErrorList.append(TUpper-TLower)
    #F0List.append(params[0])

    params, pcov=spo.curve_fit(fitLockedGamma, imageFreqMhzArr, valArr, p0=guess, bounds=bounds)
    #print(params)

    perr=np.sqrt(np.diag(pcov))
    T=1e3*calculateTemp(2.355*params[2])
    TUpper=1e3*calculateTemp(2.355*(params[2]+perr[2]))
    TLower=1e3*calculateTemp(2.355*(params[2]-perr[2]))
    TLockedList.append(T)
    TLockedErrorList.append(TUpper-TLower)
    F0List.append(params[0])

    # plt.title('SINGLE ROW OF PIXELS IN FOCUS')
    # plt.plot(imageFreqMhzArr,valArr)
    # xTest=np.linspace(imageFreqMhzArr[0],imageFreqMhzArr[-1],num=10000)
    # plt.plot(xTest,fit(xTest,*params))
    # plt.grid()
    # plt.xlabel('frequency, MHz')
    # plt.show()


F0Arr=np.asarray(F0List)
# plt.plot(zArr[yStart:yEnd],F0Arr)
# plt.grid()
# plt.show()
dVArr=3e8*F0Arr*1e6/gv.Li_D2_Freq
#plt.plot(zArr[yStart:yEnd],TList)
plt.title('Temperature vs position in focus')
plt.title('Laser jitter=2.0MHz RMS')
plt.plot(zArr[yStart:yEnd],TLockedList,label='gamma=5.87MHz')
plt.plot(zArr[yStart:yEnd],TFreeList,label='gamma free')
plt.xlabel('Vertical position in atomic beam, mm')
plt.ylabel('Temperature, mk')
plt.axvline(x=centerZ,c='r',linestyle=':',label='Atom beam center')
plt.grid()
plt.legend()
#plt.savefig('figureRun15ShutterOpen_LowJitter')
plt.show()



