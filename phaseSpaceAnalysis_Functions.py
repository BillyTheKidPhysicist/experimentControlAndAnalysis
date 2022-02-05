import warnings
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

def get_Images_Array_And_MHz_Scale_Arr(folderPath,runName):
    #Directory
    # path = "C:\Data\Runs\\8_29_21"
    os.chdir(folderPath)



    #Opening fits file and creating array. Flip the y axis such that image is oriented correctly.
    fitsFile = fits.open(runName+'.fits')
    imagesList = fitsFile[0].data
    imagesArr=imagesList.astype(float)
    if len(imagesArr.shape)!=3:
        raise Exception('Array of images is in wrong format. \n must be an array of 2d images')
    imagesArr=np.flip(imagesArr,axis=1)


    #--------------MAKE MHZ SCALE-------------------------

    #get MHS array
    infoFile=open(runName[:-4]+'Info.txt')
    imagesStartVolt=float(infoFile.readline().split(',')[1])
    imagesStopVolt=float(infoFile.readline().split(',')[1])
    #print(startVolt,stopVolt)
    dataFileName=runName[:-4]+'DAQData.csv'
    DAQData=np.loadtxt(dataFileName,delimiter=',')
    MHZMaker=MHzScale(DAQData)
    returnFitFunc=False
    MHzScaleArr=MHZMaker.make_MHZ_Scale(returnFitFunc)

    galvoVoltArr=DAQData[:, 0]
    liRefVoltArr=DAQData[:, 1]
    slope,offset=np.polyfit(galvoVoltArr, MHzScaleArr, 1) #idk why this says there is in issue in pycharm
    imageVoltageScale=np.linspace(imagesStartVolt, imagesStopVolt, num=imagesArr.shape[0])
    imagesFreqMhzArr=offset+slope*imageVoltageScale
    print('this is maybe wrong')
    return imagesArr,imagesFreqMhzArr

def subtract_Background(imagesArrOriginal,numFramesFromEnd):
    imagesArrCopy=imagesArrOriginal.copy()  #copy the array so we don't modify outside the function
    #Find background
    imageBackGround=(np.mean(imagesArrCopy[-numFramesFromEnd:], axis=0)+np.mean(imagesArrCopy[:numFramesFromEnd],
                                                             axis=0))/2  #take images from beginning and end
    #to get average of background noise
    imagesArr=imagesArrCopy-imageBackGround  #extract background.
    return imagesArr
def smooth_Images(imagesArrOriginal,filterPixelSize=1):
    imagesArrCopy=imagesArrOriginal.copy()  #copy the array so we don't modify outside the function
    #go through each image and apply a median filter to remove hot/dead pixels and cosmic rays
    for i in range(imagesArrCopy.shape[0]):
        imagesArrCopy[i]=spni.median_filter(imagesArrCopy[i], size=filterPixelSize)
    return imagesArrCopy
def find_Signal_Frequency_Bounds(imagesArr,showPlot=True,boundsFWHM_Factor=1.5):
    #automatically find the peak frame
    if imagesArr.shape[0]<15:
        warnings.warn('There may be too few images to correctly find bounds')

    xStart=50
    xEnd=-50
    yStart=50
    yEnd=-50
    trimmedImageArr=imagesArr[:, yStart:yEnd, xStart:xEnd]
    signalProfile=np.mean(np.mean(trimmedImageArr, axis=2),axis=1)
    # temp=np.mean(trimmedImageArr, axis=1)
    signalProfile=spsig.savgol_filter(signalProfile, 11, 2) #smooth the profile to assist

    guess=[np.argmax(signalProfile), np.max(signalProfile)-np.min(signalProfile), 1, 1, np.min(signalProfile)]
    x=np.arange(signalProfile.shape[0])
    params, pcov=spo.curve_fit(voigtImageFit, x, signalProfile, p0=guess)
    fwhmImage=.5346*(2*params[3])+np.sqrt(.2166*(2*params[3])**2+(params[2]*2.335)**2)

    #select the left and right image to capture most of the signal
    centerImage=params[0]
    imageLeft=int(centerImage-boundsFWHM_Factor*fwhmImage)
    imageRight=int(centerImage+boundsFWHM_Factor*fwhmImage)
    if showPlot==True:
        plt.axvline(x=imageLeft, c='black')
        plt.axvline(x=imageRight, c='black')
        plt.axvline(x=centerImage, c='green')
        plt.scatter(x, signalProfile, c='r')
        plt.plot(voigtImageFit(x, *params))
        plt.grid()
        plt.show()



# path = "C:\Data\Runs\\8_29_21"
# fileName = 'run15Far'
# get_Images_Array_And_MHz_Scale_Arr(path,fileName)