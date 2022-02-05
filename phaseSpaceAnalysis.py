import scipy.signal as spsig
import sys
from DataAnalysis import fit_Spectral_Data
import globalVariables as gv
import numpy as np
import random as ra
import os
import matplotlib.pyplot as plt
from astropy.io import fits
import scipy.special as sps
import scipy.optimize as spo
from MakeMHzScale import MHzScale
import scipy.ndimage as spni
import scipy.interpolate as spi
from UncertaintyAnalysis_Functions import find_Chi_Squared_Uncertainty
from qGaussian import qGaussianHelper
import time
import random

time_start = time.perf_counter()

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
path = "C:\Data\Runs\\9_23_21"
os.chdir(path)
fileName = 'run10Far'

#z position of the platform relative to a reference point I selected and then offset to give the correct
#distance from the combiner magnet
zpos=1000
zposOffset = 13.7
zposOffset = 0

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
# plt.axvline(x=imageLeft,c='black')
# plt.axvline(x=imageRight,c='black')
# plt.scatter(x,temp,c='r')
# plt.plot(voigtImageFit(x, *params))
# plt.grid()
# plt.show()



#Spatial Analysis Box
yStart=0
yEnd=-1
delta=-50*0
xStart=110#109-delta
xEnd=130#111-delta
image=np.mean(imagesArr[imageLeft:imageRight,yStart:yEnd,xStart:xEnd],axis=0)
# plt.imshow(image)
# plt.show()


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
# plt.suptitle('Total Focus Signal vs Frequency. Geometric broadening: '+str(lensHeating))
# plt.title('T= '+str(np.round(T,1))+' mk,   F0= '+str(np.round(freqCenter)) +' MHz')
# xTest=np.linspace(imageFreqMhzArr[0], imageFreqMhzArr[-1], num=10000)
#
# #plt.scatter(imageFreqMhzArr,valArr,c='r',label='Data',s=50.0,marker='x')
#
# plt.plot(imageFreqMhzArr,valArr,c='r')
# plt.plot(imageFreqMhzArr, spectralFit.fit_Result_Function(imageFreqMhzArr),label='Fit')
# plt.axvline(x=freqCenter,c='black',linestyle=':',label='Peak')
# plt.xlabel('frequency, MHz')
# plt.grid()
# plt.legend()
# plt.xlabel('Frequency, MHz')
# plt.ylabel('Signal Strength, AU')
# #plt.savefig('TotalTemperature_'+fileName)
# plt.show()

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
print('gamma',params[3]*2)
rMax=30 #bounds of integral in mm
rArr=np.linspace(0,rMax,num=1000)
# plt.plot(rArr,2*np.pi*rArr*voigtSpaceFit(rArr,0.0,1.0,sigma,gamma,0.0))
# plt.show()
# print("Integral",np.trapz(2*np.pi*rArr*voigtSpaceFit(rArr,0.0,1.0,sigma,gamma,0.0),x=rArr))

# plt.title('Signal vs. Transerve Position')
# plt.xlabel("Transverse position, mm")
# plt.ylabel("Signal, au")
# plt.scatter(zArr,profArr,c='r',label='Data')
# plt.plot(zArr,voigtSpaceFit(zArr,*params),label='fit')
# plt.legend()
# plt.grid()
# plt.show()











#
#
#
# #get the stacked image over the correct region
FWHM = []
FWHMSig = []
FWHMpos = []
deltaFWHM_chi = []
deltaFWHM_bootsrap = []
width=8
i = 72
offset = 16
pixel = []

while i <= (256-48-width):

    xStart=i
    xEnd=xStart+width

    yStart=0
    yEnd=-1
    temp=imagesArr[0:-1,yStart:yEnd,xStart:xEnd]

    #print('sum of images: ',np.round(np.sum(temp)/1e6,1))

    Type='mean'
    if Type=='sum':
        imageMean=np.sum(temp,axis=0) #the stacked image
    elif Type=='mean':
        imageMean=np.mean(temp,axis=0) #the stacked image
    else:
        raise Exception('INVALID')
    #print('Final image is '+Type+' of images')

    # plt.imshow(imageMean)
    # plt.title('imageMean')
    # plt.show()

#     #--------------SPATIAL PROFILE OF LASER----------------------
#     #this can be done with either a voigt fit, or a Radial Basis Function (RBF) fit. The rbf is better for data
#     #that doesn't look very voigty
#
    magnification=3 #this can be found from using a ruler, or using optics rules.
    pixelSize=magnification*4*.024 #in units of mm. binning X magnification X pixel size

    #array of the profile
    profArr=np.mean(imageMean[yStart:yEnd],axis=1) #laser profile from data

    zArr=(np.arange(0, profArr.shape[0]))*pixelSize  #each position value corresponds to the center of the pixel

    #------voigt fit
    def voigtSpaceFit(x,x0,a,sigma,gamma,b):
        v0=sps.voigt_profile(0, sigma, gamma)
        v=sps.voigt_profile(x-x0, sigma, gamma)
        return a*v/v0+b
    eps=1e-10 #small non zero number
    guess=[30.0,10000.0,1e-1,2.5,1240]
    bounds=[(-np.inf,eps,eps,eps,0.0),(np.inf,np.inf,np.inf,np.inf,np.inf)]
    params = spo.curve_fit(voigtSpaceFit, zArr, profArr, p0=guess,bounds=bounds)[0]
    fwhm=.5346*(2*params[3])+np.sqrt(.2166*(2*params[3])**2+(params[2]*2.335)**2)

    # print('voight',fwhm)
    # print('gamma',params[3]*2)

    #print('fwhm',np.round(fwhm,2),'Z0',np.round(params[0],1),'sig peak',np.round(params[1]))
    # FWHM.append(fwhm)

    FWHMSig.append(params[1])
    FWHMpos.append((zpos-(i-width/2)*0.288)/10-zposOffset)
    pixel.append(i-width/2)

    qG=qGaussianHelper()
    qGguessParams = [31.9,5000,-3.9,0.98,2.35,1.0]
    qGbounds=[(-np.inf,1e-9,-np.inf,1e-6,1e-9,1e-9),(np.inf,np.inf,np.inf,np.inf,np.inf,np.inf)]

    plt.scatter(zArr,profArr)
    plt.plot(zArr,qG(zArr,*qGguessParams))
    plt.show()

    qGparams,pcov = spo.curve_fit(qG,zArr,profArr,bounds=qGbounds,p0=qGguessParams)

    plt.scatter(zArr, profArr, c = 'r', s = 15)
    plt.plot(zArr,qG(zArr, *qGparams))
    plt.show()

    FWHMqG = qG.get_FWHM(*qGparams)
    FWHM.append(FWHMqG)
    perr = np.sqrt(np.diag(pcov))

    print('paramsOriginal',qGparams)
    print('perr',perr)

    #---qGauss FWHM Uncertainty Using Bootstrapping Method---
    if True:
        R = np.subtract(qG(zArr,*qGparams),profArr)
        Rcopy = R.copy()

        Rarr = []
        for M in range(len(Rcopy)):
            Rarr.append(Rcopy[M])

        indexArr=[]

        for iR in range(len(Rcopy)):
            indexArr.append(iR)

        sort_size=round(0.7*len(indexArr))

        sigma_i = []
        q_i = []
        eta_i = []

        for k in range(0,750):

            profArrNew = profArr.copy()
            newIndex = random.sample(indexArr,sort_size)
            newR = random.sample(Rarr,sort_size)

            for N in range(sort_size):
                profArrNew[newIndex[N]] += newR[N]

            qGparamsTemp, pcov=spo.curve_fit(qG, zArr, profArrNew, bounds=qGbounds, p0=qGguessParams)
            sigma_i.append(qGparamsTemp[3])
            q_i.append(qGparamsTemp[4])
            eta_i.append(qGparamsTemp[5])

        sigma_uncertainty = np.std(sigma_i)
        q_uncertainty = np.std(q_i)
        eta_uncertainty = np.std(eta_i)

        # plt.hist(sigma_i,bins=20)
        # plt.show()

        print('sigma uncertainty',sigma_uncertainty)
        print('q uncertainty', q_uncertainty)
        print(qGparams[3])
        FWHM_values_bootstrap = []
        for x in range(0,500):
            sigma_temp = np.random.normal(qGparams[3],sigma_uncertainty)
            q_temp = np.random.normal(qGparams[4],q_uncertainty)
            eta_temp = np.random.normal(qGparams[5],eta_uncertainty)
            FWHM_values_bootstrap.append(qG.get_FWHM(1,1,1,sigma_temp,q_temp,eta_temp))

        # plt.hist(FWHM_values,bins=20)
        # plt.show()
        FWHM_Uncertainty_bootsrap = 2*np.std(FWHM_values_bootstrap)
        print('FWHM Uncertainty BootStrapping',FWHM_Uncertainty_bootsrap)
        deltaFWHM_bootsrap.append(FWHM_Uncertainty_bootsrap)

    time_elapsed = (time.perf_counter()-time_start)
    print('time elapsed',time_elapsed)

    yvalues = qG(zArr,*qGparams)
    yvaluesshitfit = voigtSpaceFit(zArr,*params)

    # plt.plot(zArr_original,yvalues)
    # plt.plot(zArr,yvaluesshitfit,c='g')
    # plt.scatter(zArr,profArr,c='r')
    # plt.title('q Gaussian Fit')
    # plt.show()

    #--Uncertainty in Spatial FWHM From Chi Squared---

    if True:

        sigma_squared = np.square(np.subtract(yvalues,profArr))

        optimal_sigma = qGparams[3]
        optimal_q = qGparams[4]
        optimal_eta = qGparams[5]

        testSize = 0.01

        chi_sigma = []
        chi_q = []
        chi_eta = []
        test_sigma = []
        test_q = []
        test_eta = []

        for sigma_value in np.arange(optimal_sigma*(1-testSize),optimal_sigma*(1+testSize),2*testSize*optimal_sigma/150):

            def qG_fixedSigma(x,mu,a,b,q,eta):
                sigmaFixed = sigma_value
                return qG(x,mu,a,b,sigmaFixed,q,eta)

            guess_fixedSigma=[31, 500, 3.9, 1.5,1]

            # plt.plot(zArr,qG_fixedSigma(zArr,*guess_fixedSigma))
            # plt.scatter(zArr,profArr)
            # plt.show()

            bounds_fixedSigma=[(-np.inf, 1e-9, -np.inf, 1e-9,1e-9), (np.inf, np.inf, np.inf, np.inf,np.inf)]

            params_fixedSigma, pcov=spo.curve_fit(qG_fixedSigma, zArr, profArr, bounds=bounds_fixedSigma, p0=guess_fixedSigma)

            # plt.plot(zArr,qG_fixedSigma(zArr,*params_fixedSigma))
            # plt.scatter(zArr,profArr,c='r')
            # plt.title('Fit for Test Parameters')
            # plt.show()

            R_sigma = np.subtract(profArr,qG_fixedSigma(zArr,*params_fixedSigma))
            R2_sigma = np.square(R_sigma)

            chi_squared_sigma=np.sum(np.divide(R2_sigma,sigma_squared))

            test_sigma.append(sigma_value)
            chi_sigma.append(chi_squared_sigma)

        for q_value in np.arange(optimal_q*(1-testSize),optimal_q*(1+testSize),2*testSize*optimal_q/150):

            def qG_fixedq(x,mu,a,b,sigma,eta):
                qFixed = q_value
                return qG(x,mu,a,b,sigma,qFixed,eta)

            guess_fixedq=[31.9, 800, 3.9, 2.5,1]

            bounds_fixedq=[(-np.inf, 1e-9, -np.inf, 1e-9,0), (np.inf, np.inf, np.inf, np.inf,np.inf)]

            params_fixedq, pcov=spo.curve_fit(qG_fixedq, zArr, profArr, bounds=bounds_fixedq, p0=guess_fixedq)

            # plt.plot(zArr,qG_fixedSigma(zArr,*params_fixedSigma))
            # plt.scatter(zArr,profArr)
            # plt.show()

            R_q = np.subtract(profArr,qG_fixedq(zArr,*params_fixedq))
            R2_q = np.square(R_q)
            chi_squared_q = np.sum(np.divide(R2_q,sigma_squared))

            test_q.append(q_value)
            chi_q.append(chi_squared_q)

        for eta_value in np.arange(optimal_eta*(1-20*testSize),optimal_eta*(1+15*testSize),2*testSize*optimal_eta/150):

            def qG_fixedeta(x,mu,a,b,sigma,q):
                etaFixed = eta_value
                return qG(x,mu,a,b,sigma,q,etaFixed)

            guess_fixedeta=[31.9, 800, 3.9, 2.5,2]

            bounds_fixedeta=[(-np.inf, 1e-9, -np.inf, 1e-9,0), (np.inf, np.inf, np.inf, np.inf,np.inf)]

            params_fixedeta, pcov=spo.curve_fit(qG_fixedeta, zArr, profArr, bounds=bounds_fixedeta, p0=guess_fixedeta)

            # plt.plot(zArr,qG_fixedeta(zArr,*params_fixedeta))
            # plt.scatter(zArr,profArr)
            # plt.show()

            R_eta = np.subtract(profArr,qG_fixedeta(zArr,*params_fixedeta))
            R2_eta = np.square(R_eta)
            chi_squared_eta = np.sum(np.divide(R2_eta,sigma_squared))

            test_eta.append(eta_value)
            chi_eta.append(chi_squared_eta)

        sigma_uncertainty2=find_Chi_Squared_Uncertainty(test_sigma, chi_sigma, False)
        q_uncertainty2 = find_Chi_Squared_Uncertainty(test_q, chi_q,False)
        eta_uncertainty2 = find_Chi_Squared_Uncertainty(test_eta,chi_eta,False)

        FWHM_values=[]
        for x in range(0, 500):
            sigma_temp=np.random.normal(qGparams[3],sigma_uncertainty2)
            q_temp=np.random.normal(qGparams[4],q_uncertainty2)
            eta_temp = np.random.normal(qGparams[5],eta_uncertainty2)
            FWHM_values.append(qG.get_FWHM(1, 1, 1, sigma_temp, q_temp,eta_temp))

        # plt.hist(FWHM_values,bins=25)
        # plt.show()
        FWHM_Uncertainty_Chi= 2*np.std(FWHM_values)
        print('FWHM Uncertainty Chi Squared', FWHM_Uncertainty_Chi)

        deltaFWHM_chi.append(FWHM_Uncertainty_Chi)

    i = i+width + offset

FWHM_array = [FWHMpos, FWHM, deltaFWHM_chi,deltaFWHM_bootsrap]
#
# txtfile = open("FWHM1.txt", "r")
# np.savetxt("FWHM1.txt",FWHM_array)

print(FWHMpos)
print(FWHM)
print('deltaFWHM Chi',deltaFWHM_chi)
print('deltaFWHM Boot',deltaFWHM_bootsrap)

fig, ax = plt.subplots()
ax.errorbar(FWHMpos,FWHM,yerr=deltaFWHM_chi,fmt='.k',marker='o',mfc='w',color='blue')
plt.show()

fig, ax = plt.subplots()
ax.errorbar(FWHMpos,FWHM,yerr=deltaFWHM_bootsrap,fmt='.k',marker='o',mfc='w',color='blue')
plt.show()

# plt.scatter(FWHMpos,FWHM)
# plt.title('FWHM of Voigt vs Position for '+str(fileName))
# plt.xlabel("Distance from Face of Output of Magnet (mm)")
# plt.ylabel("FWHM of Voigt (mm)")
# #plt.savefig('FWHM_'+fileName)
# plt.grid()
# plt.show()


plt.scatter(pixel,FWHM)
plt.title('FWHM of Voigt vs Position for '+str(fileName))
plt.xlabel("Pixel #")
plt.ylabel("FWHM (mm)")
plt.grid()
#plt.savefig('FWHMSignal_'+fileName)
plt.show()

