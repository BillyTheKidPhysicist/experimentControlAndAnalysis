import scipy.signal as spsig
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

#analyze the big fat wide focus. Also make a pretty picture of it if so desired


#------voigt fit
def voigtSpaceFit(x,x0,a,sigma,gamma,b):
    v0=sps.voigt_profile(0, sigma, gamma)
    v=sps.voigt_profile(x-x0, sigma, gamma)
    return a*v/v0+b
def voigtImageFit(x, x0, a, sigma, gamma, b):
    v0=sps.voigt_profile(0, sigma, gamma)
    v=sps.voigt_profile(x-x0, sigma, gamma)
    return a*v/v0+b


#Directory
path = "C:\Data\Runs\\6_23_21"
os.chdir(path)
fileName = 'run18Far'
#Opening fits file and creating array. Flip the y axis such that image is oriented correctly.
fitsFile = fits.open(fileName+'.fits')
imagesList = fitsFile[0].data 
imagesArr=imagesList.astype(float)
imagesArr=np.flip(imagesArr,axis=1)
trimVal=3e3
imagesArr[imagesArr>trimVal]=trimVal


for i in range(imagesArr.shape[0]):
    imagesArr[i]=spni.gaussian_filter(imagesArr[i], 0.5)


#crop the image to the focus to remove noise
xStart=0
xEnd=-1
yStart=0
yEnd=800
temp=imagesArr[:,yStart:yEnd,xStart:xEnd]


temp=np.mean(temp,axis=2)
temp=np.mean(temp,axis=1)




guess=[np.argmax(temp),np.max(temp)-np.min(temp),1,1,np.min(temp)]
x=np.arange(temp.shape[0])
params,pcov=spo.curve_fit(voigtImageFit, x, temp, p0=guess)
fwhmImage=.5346*(2*params[3])+np.sqrt(.2166*(2*params[3])**2+(params[2]*2.335)**2)

#select the left and right image to capture most of the signal
centerImage=params[0]
imageLeft=int(centerImage-2*fwhmImage)
imageRight=int(centerImage+2*fwhmImage)

plt.axvline(x=imageLeft,c='black')
plt.axvline(x=imageRight,c='black')
plt.scatter(x,temp,c='r')
plt.plot(voigtImageFit(x, *params))
plt.grid()
plt.show()



#this start and end value captures most of the beam in frequency space. It is used to make the single image that is the mean
#of a stack of iamges between imStart:imEnd
offset=0
imStart=imageLeft
imEnd=imageRight
imageMean=np.mean(imagesArr[imStart:imEnd],axis=0)






#visualize the results
binning=6
magnification=2.77
print('BINNING IS SET TO '+str(binning))
pixelSize=magnification*binning*.025 #in mm
print(pixelSize)

#-----------prepare the mean image-----------------------------------------
#Position in mm. Zero corresponds to leftmost side of the laser beam. The start of the viewing region is at 8mm
powerPosArr = np.asarray([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85])-8 #subtract 8 to move from laser
#frame to image frame
#Measuring the power with the power meter with no small aperature over the sensor.
powerArr = np.asarray([37.2, 33.2, 32.3, 33, 34.1, 35.9, 38.2, 38.3, 37, 35, 33, 31.6, 30, 27, 26.3, 25, 23.2, 22.3])
powerArr=powerArr/np.max(powerArr) #normalize data
powerCorrectionFunc=spi.Rbf(powerPosArr, 1/powerArr, smooth=.1)
# test=np.linspace(0,85)
# plt.plot(test,powerCorrectionFunc(test))
# plt.scatter(powerPosArr,1/powerArr)
# plt.show()

def adjust_Image_To_Laser_Power(image):
    #adjust each column in image by a factor which is equal to the inverse of the laser power (normalized to 1)
    #at that location. The goal is to try and remove the effect of laser power. This is done by making a correction
    #matrix with the same dimensions as the image and multiplyuing the image with it. A column in the correction matrix
    #has the same correction factor in it up and down
    powerCorrection=powerCorrectionFunc(np.arange(image.shape[1])*pixelSize)
    powerCorrectionMatrix=np.tile(powerCorrection,[image.shape[0],1])
    return image*powerCorrectionMatrix
def adjust_Image_To_Peak_Brightness_Each_Column(image):
    peakSignalArr=np.max(image,axis=0)
    print(np.argmax(image,axis=0))
    peakSignalArr=peakSignalArr/peakSignalArr.max()
    correctionArr=1/peakSignalArr
    # correctionArr=spsig.savgol_filter(correctionArr,2*int(correctionArr.shape[0]*.25)//2+1,1) #if you want to smooth it
    #. You can play with the parameters, it's rather qualitative
    pixelArr=np.arange(image.shape[1])
    # plt.plot(pixelArr,correctionArr)
    # plt.plot(pixelArr,1/peakSignalArr)
    # plt.show()
    correctionMatrix=np.tile(correctionArr,[image.shape[0],1])
    return image*correctionMatrix


imageBackGround=(np.mean(imagesArr[-5:],axis=0)+np.mean(imagesArr[:5],axis=0))/2 #take images from beginning and end
#to get average of background noise
imageBackGround=imageBackGround[yStart:yEnd,xStart:xEnd]
#






prettyImage=imageMean[yStart:yEnd,xStart:xEnd].copy() #crop image to only the focus to remove noise. copy to not mess
#with the original image by accident if using the whole image
prettyImage=prettyImage-imageBackGround
# prettyImage=spni.gaussian_filter(prettyImage, 1.0)
# prettyImage=adjust_Image_To_Laser_Power(prettyImage)
prettyImage=adjust_Image_To_Peak_Brightness_Each_Column(prettyImage)
cutoff=30
prettyImage[prettyImage<cutoff]=cutoff
cutoff2=180
prettyImage[prettyImage>cutoff2]=cutoff2
plt.imshow(prettyImage)#,cmap='gray')
plt.show()







#now go through the focus and get the transvese profile

#first define function to quantify the width
def width_Function(y):
    x=np.arange(y.shape[0])*pixelSize
    x0Guess=x[np.argmax(y)]
    aGuess=np.max(y)-np.min(y)
    sigmaGuess=1.0
    gammaGuess=1.0
    bGuess=x[np.argmin(y)]
    guess=[x0Guess,aGuess,sigmaGuess,gammaGuess,bGuess]
    eps=1e-10
    bounds=[(-np.inf, eps, eps, eps, 0.0), (np.inf, np.inf, np.inf, np.inf, np.inf)]
    try:
        params,pcov=spo.curve_fit(voigtSpaceFit,x,y,p0=guess,bounds=bounds)
    except:
        print('failed to fit')
        return np.nan
    fwhm=.5346*(2*params[3])+np.sqrt(.2166*(2*params[3])**2+(params[2]*2.335)**2)
    # print(params[2:4],fwhm)
    # plt.scatter(x,y,c='r')
    # plt.plot(x,voigtSpaceFit(x,*params))
    # plt.show()
    return fwhm,params[1]



yUpper=yStart #beginning of transverse plot in the vertical direction
yLower=yEnd #ending of transverse plot in the vertical direction

numColumns=prettyImage.shape[1]-1
print('numcolumns',numColumns)
columnSteps=6 #aggregate this many columns into one
boxSize = np.round(columnSteps*pixelSize,2)
widthList=[]
xWidth = numColumns*pixelSize
sigList=[]
if numColumns%columnSteps!=0:
    print('column clumping amount should be a divisor of number of columns with no remainder')
numColumns=(numColumns//columnSteps)*columnSteps
for i in range(0,numColumns,columnSteps):
    y=np.sum(prettyImage[:,i:i+columnSteps],axis=1)
    width,sig=width_Function(y)
    widthList.append(width)
    sigList.append(sig)
totalSignal=np.mean(prettyImage)
xArr = np.linspace(0,xWidth,num=int(numColumns/columnSteps))
# zArr=np.arange(0,numColumns/columnSteps)*pixelSize

plt.scatter(xArr,sigList)
plt.xlabel('Distance from leftward edge of frame (cm)')
plt.ylabel('Signal Strength of Voigt for Each Box')
plt.title(fileName+' ,Box Size in x of %s mm' %boxSize)
plt.show()

plt.scatter(xArr,widthList)
plt.title(fileName+' ,Box Size in x of %s mm' %boxSize)
plt.xlabel('Distance from leftward edge of frame (cm)')
plt.ylabel('FWHM for each Box (mm)')
# plt.savefig(fileName)
plt.show()


