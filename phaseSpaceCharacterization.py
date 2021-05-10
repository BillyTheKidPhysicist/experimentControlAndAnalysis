
import globalVariables as gv
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import fits
import scipy.signal as spsig
import scipy.special as sps
import scipy.optimize as spo
from MakeMHzScale import MHzScale
import scipy.ndimage as spni
from DataAnalysis import fitPixelData, spectralProfile
import scipy.interpolate as spi


#----------------IMPORT AND RANGLE IMAGES INTO RIGHT SHAPE-----------
#Directory
path = "C:\Data\Runs\\3_22_21"
os.chdir(path)
fileName = 'run3Far'


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
for i in range(imagesArr.shape[0]):
    imagesArr[i]=spni.gaussian_filter(imagesArr[i], .5) #slightly smooth each image


#this start and end value captures 90% of the beam in frequency space. It is used to make the single image that is the mean
#of a stack of iamges between imStart:imEnd
offset=-5
imStart=70+offset
imEnd=100+offset
imageMean=np.mean(imagesArr[imStart:imEnd],axis=0)







#crop the image to the focus to remove noise
xStart=81
xEnd=82
yStart=77
yEnd=98




#visualize the results


# #plot what the mean image looks like
temp=imageMean[yStart:yEnd,xStart:xEnd] #crop image to only the focus to remove noise
plt.imshow(temp)
plt.show()

# #plot the mean of the focus over the restricted images
temp=imagesArr[imStart:imEnd,yStart:yEnd,xStart:xEnd]
temp=np.mean(temp,axis=2)
temp=np.mean(temp,axis=1)
plt.scatter(np.arange(temp.shape[0]),temp)
plt.plot(np.arange(temp.shape[0]),temp)
plt.grid()
plt.show()



#--------------SPATIAL PROFILE OF LASER----------------------
#this can be done with either a voigt fit, or a Radial Basis Function (RBF) fit. The rbf is better for data
#that doesn't look very voigty

#expand the image region vertically
ya=yStart-23
yb=yEnd+23

magnification=2.77 #this can be found from using a ruler, or using optics rules.
pixelSize=magnification*5*.025 #in units of mm. binning X magnification X pixel size

#array of the profile
profArr=np.mean(imageMean[ya:yb,xStart:xEnd],axis=1) #laser profile from data


zArr=(np.arange(0,imagesArr[0].shape[0]))*pixelSize #each position value corresponds to the center of the pixel
print((zArr.max()-zArr.min())/2)



#two fit options, voigt or RBF

#------voigt fit
def voigtSpaceFit(x,x0,a,sigma,gamma,b):
    v0=sps.voigt_profile(0, sigma, gamma)
    v=sps.voigt_profile(x-x0, sigma, gamma)
    return a*v/v0+b
eps=1e-10 #small non zero number
guess=[30.0,290.0,1e-1,2.5,1240]
bounds=[(-np.inf,eps,eps,eps,0.0),(np.inf,np.inf,np.inf,np.inf,np.inf)]
params, pcovLoopGammaFree=spo.curve_fit(voigtSpaceFit, zArr[ya:yb], profArr, p0=guess,bounds=bounds,)
fwhm=.5346*(2*params[3])+np.sqrt(.2166*(2*params[3])**2+(params[2]*2.335)**2)
print('sigma',params[2],'gamma',params[3], 'mm')
print('fwhm',np.round(fwhm,3),'Z0',np.round(params[0],1),'sig peak',np.round(params[1]))


zArrFine=np.linspace(zArr[ya:yb].min(),zArr[ya:yb].max(),num=1000) #very fine z array
voigtFitArrFine=voigtSpaceFit(zArrFine,*params) #very fine fit array


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
#plt.savefig('Spatial_'+fileName)
plt.show()
plt.close('all')



#-----------------fit the veocity dependence on space----------------------
#--------------------------------------------------------
def spectralProfileFitter(x, x0, a, sigma, gamma, b, singlePeak=False):
    #x: frequency, MHz
    #x0: center frequency, MHz
    #a: peak height, au
    #sigma: standard deviation of gaussian contribution, MHz
    #gamma: HWHM of lorentzian, MHz
    #singlePeak: If true, assume optical pumping, if False use multipeak fit
    laserJitter=2/2.355
    sigma=np.sqrt(sigma**2+laserJitter**2)
    if singlePeak==True:
        v0=sps.voigt_profile(0,sigma,gamma)
        v=sps.voigt_profile(x-x0,sigma,gamma)
        return a*v/v0 +b
    else:
        val=spectralProfile(x,a,x0,b,dF0=sigma*2.355,gamma=2*gamma) #this function uses FWHM for linewidhts and constants are
        #ordered differently
        return val



zVals=[]
sigmaVals=[]
F0Vals=[]
#-----------------get the peak nice and centered-----------------------
yStart-=0
yEnd+=-2
profArr=np.mean(imageMean[yStart:yEnd,xStart:xEnd],axis=1) #laser profile from data
plt.plot(zArr[yStart:yEnd],profArr)
plt.axvline(x=centerZVoigt+params[3])
plt.axvline(x=centerZVoigt+2*params[3])
plt.axvline(x=centerZVoigt-params[3])
plt.axvline(x=centerZVoigt-2*params[3])
plt.grid()
plt.show()



#-------------fit the whole thing with a voigt to frequency----------------------
imageMeanArr=imagesArr[:,yStart:yEnd,xStart:xEnd]
imageMeanArr=np.mean(imageMeanArr,axis=2)
imageMeanArr=np.mean(imageMeanArr,axis=1)
guess=[0, 100.0, 4.0, 5, 1200.0]
bounds=((-np.inf, 0, 0, 5.87/2.0, 0), (np.inf, np.inf, np.inf, np.inf, np.inf))
params, pcov=spo.curve_fit(spectralProfileFitter, imageFreqMhzArr, imageMeanArr, p0=guess, bounds=bounds)
print(params)
f0=(3e8/671e-9)/1e6
dfTodv=3e8/(f0)#frequency shift to velocity shiftb
fwhm=.5346*(2*params[3])+np.sqrt(.2166*(2*params[3])**2+(params[2]*2.335)**2)
print('fwhm',fwhm*dfTodv, 'sigma', params[2]*dfTodv,'m/s')
xPlot=np.linspace(imageFreqMhzArr[0],imageFreqMhzArr[-1],num=10000)


plt.scatter(imageFreqMhzArr,imageMeanArr,alpha=.5)
plt.plot(imageFreqMhzArr,imageMeanArr,alpha=.5)
plt.plot(xPlot,spectralProfileFitter(xPlot,*params))
plt.grid()
plt.show()








#--------------fit each pixel with a voigt-------------------------
for i in range(yStart,yEnd):
    #print(i)
    valArr=np.mean(imagesArr[:, i, xStart-1:xEnd+1], axis=1)
    guess=[0, 100.0, 4.0, 5, 1200.0]
    bounds=((-np.inf, 0, 0, 5.87/2.0, 0), (np.inf, np.inf, np.inf, np.inf, np.inf))
    params,pcov=spo.curve_fit(spectralProfileFitter,imageFreqMhzArr,valArr,p0=guess,bounds=bounds)
    # plt.plot(imageFreqMhzArr,valArr)
    # plt.plot(imageFreqMhzArr,spectralProfileFitter(imageFreqMhzArr,*params))
    # plt.show()
    sigmaVals.append(params[2])
    zVals.append(zArr[i])
    F0Vals.append(params[0])
    #plt.close('all')
    #plt.plot(valArr)
    #plt.show()
f0=(3e8/671e-9)/1e6
dfTodv=3e8/(f0)#frequency shift to velocity shift
sigmaVals=np.asarray(sigmaVals)*dfTodv
F0Vals=np.asarray(F0Vals)*dfTodv


#slope of left half is ~-.50 (m/s)/mm. Minimum at focus is about 4m/s
plt.plot(zVals,sigmaVals*2.355)
plt.grid()
plt.show()


m,b=np.polyfit(zVals,F0Vals,1)
print('slope is (m/s)/mm ',m) #slope is (m/s)/mm  -0.5395832797041357

plt.plot(zVals,F0Vals)
plt.grid()
plt.show()
