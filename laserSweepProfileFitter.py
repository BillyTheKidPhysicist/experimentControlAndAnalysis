import matplotlib.pyplot as plt
import numpy as np
from DAQClass import DAQPin
from CameraClass import Camera
import globalVariables as gv
from tqdm import tqdm
import time
import sys
import scipy.special as sps
import scipy.optimize as spo
def voigt(x,x0,a,b,sigma,gamma):
  v0=sps.voigt_profile(0,sigma,gamma)
  return b+a*sps.voigt_profile(x-x0,sigma,gamma)/v0





minVolt=-.9
maxVolt=.1
numPoints=25

#x0,y0,z0 and only used to make the file name for saving. They are the position of the stage
x0=110
z0=95
y0=90
save=False #wether to save files or not. Set to False for not saving files and to only show the plots


camera=Camera('FAR',1000,binx=16,biny=2,imageParams=[0,1024,0,1024],temp=-25)

voltArr=np.linspace(minVolt,maxVolt,num=numPoints)

gv.begin_Sound()
galvoOut=DAQPin(gv.galvoOutPin)
#galvoOut.write(minVolt)
shape=camera.aquire_Image().shape #this serves 2 purposes. First, with Ximea and FLI there seems to be an advantage
#to taking a fake image first. Often the first image is different than the following images. I have no idea why. Second,
#this lets me get the shape of the array. You can't get that from what the user set here because Camera may change it
#to respect binning and image boundaries.
numDarkImages=3 #number of dark images to take. I set to no less than 3
imgDark=np.zeros(shape) #do not make the mistake of setting this value to 0. It caused alot of issues. There is some
#subtelty going on here. This method works as expected
for i in range(numDarkImages):
    img=camera.aquire_Image()
    imgDark+=img
imgDark=imgDark/numDarkImages #need to average them!


imageStacked=np.ones(shape)
for volt in tqdm(voltArr):
    #galvoOut.write(volt)
    img=camera.aquire_Image()
    img=img-imgDark #subtract the dark image
    imageStacked+=img
camera.close()
galvoOut.close()
gv.finished_Sound()
name='z'+str(z0)+'y'+str(y0)+'x'+str(x0)
if save==True:
    np.savetxt('npImageFile'+name,imageStacked)
    plt.savefig(name+'Img')
y=np.sum(imageStacked,axis=1)
x=np.linspace(0,y.shape[0]-1,num=y.shape[0])

#plt.plot(x,y)
#plt.show()

#params =x0,a,b,sigma,gamma
x0Guess=np.argmax(y)
a0Guess=np.max(y)-np.min(y)
b0Guess=np.min(y)
widthGuess=y.shape[0]/10
p0=[x0Guess,a0Guess,b0Guess,widthGuess,widthGuess]
print(p0)
bounds=([0,0,-1000,0,0],[250,2000,1000,100,100])
params=spo.curve_fit(voigt,x,y,p0=p0)[0]
sigma=params[3]
gamma=params[4]

#find the FWHM. I can use two methods here.
#method 1
#fg=2*sigma*np.sqrt(2*np.log(2))
#fl=2*gamma
#fwhm=.5346*fl+np.sqrt(.2166*fl**2+fg**2)

#method 2

xTemp=np.linspace(x[0],x[-1],num=10000)
yTemp=voigt(xTemp,*params)
yHalf=(np.max(yTemp)-np.min(yTemp))/2
HWHM=x[np.argmax(yHalf>yTemp)]
FWHM=2*HWHM




plt.close('all')
plt.plot(x,y,label='data')
plt.plot(x,voigt(x,*params),label='voigtFit')
plt.title('fwhm: '+str(FWHM))
plt.grid()
plt.xlabel('y dimension, pixels')
plt.legend()
if save==True:
    plt.savefig(name)
plt.show()
