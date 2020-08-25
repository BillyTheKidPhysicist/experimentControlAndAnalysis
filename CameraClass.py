import sys
import globalVariables as gv
import matplotlib.pyplot as plt
import pyfli as fli
import time
import progressbar
import numpy as np
from ximea import xiapi as xi
import threading
class Camera:
    def __init__(self,camName,expTime,imageParams=None,bin=1,temp=-20):
        #camName: Name of camera, or rather position. Far field vs near field
        #expTime: exposure time of camera in ms
        #imageParam: A list of image paramters, [x1,x2,y1,y2], where x1 is start of image region along x, x2 is end of
                #image region along x. y1 is start of image region along y, y2 is end of image region along y.
                #if none then use the entire sensor
        self.camName=camName
        self.expTime=expTime #convert from milliseconds to microseconds
        self.imageParams=imageParams
        self.bin=bin
        self.temp=temp
        self.cameraIndex=None #used for FLI camera. For some annoying reason FLI libraries use indices to track the camera
                #instead of a camera object
        self.cameraObject=None #ximea API does the more reasonable thing of using a python object so as of now this is
                #only for the ximea camera.
        self.initialize_Camera()
    def initialize_Camera(self):
        if self.camName=='FAR':
            self._initialize_Camera_FLI()
        elif self.camName=='NEAR':
            self._initialize_Camera_Ximea()
        else:
            raise Exception('NO VALID CAMERA NAME PROVIDED')
    def _initialize_Camera_Ximea(self):
        self.cameraObject=xi.Camera()
        self.cameraObject.open_device()
        self.cameraObject.set_imgdataformat('XI_MONO16')  #return data with the greatest depth possible. Otherwise data
        #is forced to fit into an 8 bit array
        self.cameraObject.set_exposure(int(self.expTime*1E3)) #exposure time is in microseconds for ximea camera

        #turn off every feature that I can. Who knows whats happening if these are on? I'm not even confident in this
        #camera with these off
        self.cameraObject.set_gain(0)
        self.cameraObject.disable_bpc()
        self.cameraObject.disable_ffc()
        self.cameraObject.disable_aeag()
        self.cameraObject.disable_LUTEnable()
        #----disable obnoxious image clamping thing----
        self.cameraObject.set_param('device_unit_register_selector', 0x0222)
        self.cameraObject.set_param('device_unit_register_value', 0xF0)

        self._catch_And_Fix_Errors()
        x1,x2,y1,y2=self.imageParams
        xWidth=x2-x1
        yWidth=y2-y1
        yMax=self.cameraObject.get_height_maximum()
        self.cameraObject.set_width(xWidth)
        self.cameraObject.set_height(yWidth)
        self.cameraObject.set_offsetX(x1)
        self.cameraObject.set_offsetY(yMax-yWidth-y1) #image coordinates are top left, but I want from bottom left
        self.cameraObject.start_acquisition()

    def _catch_And_Fix_Errors(self):
        x1, x2, y1, y2=self.imageParams
        if self.camName=='FAR':
            if self.bin<0 or self.bin>16:
                print('----------ERROR-------------')
                print('VALID RANGE OF BINNING FOR FLI CAMERA IS 1 TO 16')
                gv.error_Sound()
                sys.exit()
            if self.expTime<50:
                print('----------ERROR-------------')
                print('FLI CAMERA EXPOSURE MUST BE AT LEAST 50 ms')
                gv.error_Sound()
                sys.exit()
        #now adress issues common to both cameras
        if self.camName=='NEAR':
            xMax=self.cameraObject.get_width_maximum()
            yMax=self.cameraObject.get_height_maximum()
        elif self.camName=='FAR':
            xMax=1024
            yMax=1024
        else:
            raise Exception('NO VALID CAMERA NAME PROVIDED')

        if x1<0 or x2<0 or y1<0 or y2<0:
            print('NO IMAGE DIMENSION VALUES CAN BE NEGATIVE')
            gv.error_Sound()
            sys.exit()
        if x2>xMax or y2>yMax:
            print('-------WARNING-----------')
            print('ATTEMPTING TO USE AN IMAGE REGION THAT IS LARGER THAN AVAILABLE SIZE \n'
                  'MAXIMUM DIMENSIONS ARE: '+str(xMax)+','+str(yMax))
            print('CONFLICTING DIMENSION WILL BE CHANGED TO MAXIMUM')
            print('-----------END OF WARNING------------')
            gv.warning_Sound()
            if x2>xMax:
                x2=xMax
            if y2>yMax:
                y2=yMax
        if x2<=x1 or y2<=y1:
            print('END OF IMAGE REGION CANNOT BE BEFORE OR AT BEGINNING BEGINNING')
            gv.error_Sound()
            sys.exit()
        if x1==x2 or y1==y2:
            raise Exception('IMAGE BOUNDS CANNOT BE THE SAME. X1/=X2 AND Y1/=Y2')


        #make image size and binning line up correctly
        if self.camName=='NEAR':
            x1New=self.bin*8*(x1//(self.bin*8))
            x2New=self.bin*8*(x2//(self.bin*8))
            y1New=self.bin*8*(y1//(self.bin*8))
            y2New=self.bin*8*(y2//(self.bin*8))
            if x1New!=x1 or x2New!=x2 or y1New!=y1 or y2New!=y2:
                print('----------WARNING-----------')
                print('Image region for Ximea needs to be both even multiple of 8 as well as an even')
                print('multiple of the bin size. Adusting now.')
                print('Previous value [x1,x2,y1,y2]: ', [x1, x2, y1, y2])
                x1=x1New
                x2=x2New
                y1=y1New
                y2=y2New
                print('New values are now [x1,x2,y1,y2]: ', [x1, x2, y1, y2])
                print('---------------END OF WARNING---------------')
                gv.warning_Sound()
        elif self.camName=='FAR':
            x1New=self.bin*(x1//self.bin)
            x2New=self.bin*(x2//self.bin)
            y1New=self.bin*(y1//self.bin)
            y2New=self.bin*(y2//self.bin)
            if x1New!=x1 or x2New!=x2 or y1New!=y1 or y2New!=y2:
                print('----------WARNING-----------')
                print('Image region needs to be an even multiple of bin size. Adjusting now.')
                print('Previous value [x1,x2,y1,y2]: ', [x1, x2, y1, y2])
                x1=x1New
                x2=x2New
                y1=y1New
                y2=y2New
                print('New values are now [x1,x2,y1,y2]: ', [x1, x2, y1, y2])
                print('---------------END OF WARNING---------------')
                gv.warning_Sound()
        #after massaging the parameters to work with the binning and other restrictions, it's possible that all value
        #could now equal each other
        if x1==x2 or y1==y2:
            print('------------ERROR-------------')
            print('AFTER ADJUSTING IMAGE PARAMETERS TO FIT CONSTRAINTS X OR Y IMAGE DIMENSIONS NOW EQUAL EACH OTHER')
            print('PICK DIFFERENT IMAGE PARAMETERS THAT ARE SEPERATED BY MORE THAN A MULTIPLE OF THE PRODUCT OF ')
            print('BIN SIZE AND MINIMUM SEPERATION')
            gv.error_Sound()
            sys.exit()
        self.imageParams=[x1,x2,y1,y2] #change image params list to corrected values

    def _initialize_Camera_FLI(self):
        # I don't set the image region because it's a little complicatd with the binning. Instead I just take a full image
        #and format before returning




        #cameraModel=fli.FLIList('usb', 'camera')#[0][1]
        try:
            self.cameraIndex=fli.FLIOpen('flipro0', 'usb', 'camera')
        except:
            raise Exception('NO FLI CAMERA CONNECTED, OR INCORRECT FLI CAMERA CONNECTED. READ COMMENTS BELOW IN CODE FOR MORE INFO')
            #For some reason when trying to get a list of connected FLI cameras and their names if there are none
            #connected then a silent error occurs and python exits. This does not occur if I do FLIOpen first, but it
            #depends on me knowing the name of the camera which I think is always 'flipro0' for the first one. This error
            #will be thrown if that isn't true as well, such as using a different FLI camera, or some other change
            #in ports or something. Wish it didn't have to be like this but this library is a home made wrapper
        fli.setTemperature(self.cameraIndex, self.temp)
        currentTemp=fli.readTemperature(self.cameraIndex, 'internal')
        print('Opened FLI Camera . Current temperature is '+str(currentTemp)+' Celcius.')
        if currentTemp>self.temp+.5:
            print('Current temperature is above target temperatute of '+str(self.temp)+' Celcius. Initiating cooling.')
            print('It takes about 5 minutes to cool from 20 celcius to -20 celcius')
            fli.setTemperature(self.cameraIndex, self.temp)
            loop=True
            while loop==True:
                currentTemp=fli.readTemperature(self.cameraIndex, 'internal')
                time.sleep(5)
                if currentTemp>self.temp:
                    print('Current temperature: '+str(currentTemp)+' Celcius')
                if currentTemp<=self.temp:
                    print('Target temperature reached')
                    loop=False
        fli.setExposureTime(self.cameraIndex, self.expTime)



        self._catch_And_Fix_Errors()
        fli.setVBin(self.cameraIndex,self.bin)
        fli.setHBin(self.cameraIndex,self.bin)
    def aquire_Image(self):
        if self.camName=='FAR':
            #the image needs to be massaged to display how I want it to. For an image with dimension (M,N) M defined rows
            #and N defines columns. if I want to select a sub image from this I need to not confuse M with x and N with y
            #because it is the opposite actually
            img=self._aquire_Image_FLI()
            img= img[:int(img.shape[0]/self.bin),:int(img.shape[1]/self.bin)] #if the image is binned then the FLI camera
                    #returns the entire 1024x1024 image, but with a section of it containing the actual image.
            img=np.rot90(img) #image is sideways relative to mounting screw, needs to be rotated
            x1, x2, y1, y2=self.imageParams
            img=img[(1024-y2)//self.bin:(1024-y1)//self.bin,x1//self.bin:x2//self.bin]
            return img
        elif self.camName=='NEAR':
            img=self._aquire_Image_Ximea()
            return img
        else:
            raise Exception('NO VALID CAMERA NAME PROVIDED')

    def _aquire_Image_Ximea(self):
        imgObject=xi.Image()
        self.cameraObject.get_image(imgObject)
        img=imgObject.get_image_data_numpy()
        if self.bin!=1:
            img=self._bin_Image(img,self.bin)
        return img
    def _aquire_Image_FLI(self):

        fli.exposeFrame(self.cameraIndex)
        while fli.getExposureStatus(self.cameraIndex)>0: #when there is 0 milliseconds remaining the exposure is done
            time.sleep(.01) #I don't think I want to query the camera for exposure status millions of time a second
        img=fli.grabFrame(self.cameraIndex)

        return img

    def _bin_Image(self,image, binSize):
        #ximea images have to be manually binned unfortunately...
        m, n=image.shape
        return image.reshape(m//binSize, binSize, n//binSize, binSize).sum(3).sum(1)
    def close(self):
        if self.camName=='NEAR':
            #self.cameraObject.stop_acquisition()
            self.cameraObject.close_device()
#currentCam=xi.Camera()
#currentCam.open_device()
#print(currentCam.get_device_name())
#currentCam.close_device()
#
#import time
#import matplotlib.pyplot as plt
#print(fli.FLIList('usb','camera'))
#cam = fli.FLIOpen('flipro0', 'usb', 'camera')
#print(fli.readTemperature(cam,'internal'))
#fli.setTemperature(cam,0)
##print(fli.readTemperature(cam,'internal'))
#
#expTime=.1
#print(expTime)
#fli.setExposureTime(cam,expTime*1E3)
#fli.exposeFrame(cam)
#time.sleep(1)
#
#img=fli.grabFrame(cam)
#print(img.mean())
#print(img.sum())
##fli.endExposure(cam)
#print(type(img))
#plt.imshow(img)
#plt.show()
#import scipy.interpolate as spi
#import scipy.signal as sps
import astropy.io as astroio
#camera=Camera('FAR',250,imageParams=[0,500,0,500],bin=1)
#img=camera.aquire_Image()
#plt.imshow(img)
#plt.show()

from astropy.io import fits
#camera.close()

#hdu=fits.PrimaryHDU(img)
#hdul=fits.HDUList([hdu])
#hdul.writeto('test.fits')
#
#
#
#yFWHMList=[]
#for i in range(10):
#    img=camera.aquire_Image()
#    y=img[490,430:465]
#    y=y/y.max()
#
#
#    x=np.arange(0,y.shape[0])
#    ySmooth=sps.savgol_filter(y,7,2)
#    yFit=spi.Rbf(x,ySmooth,smooth=.1)
#    xArr=np.linspace(0,np.max(x),num=10000)
#
#
#    yFit=yFit(xArr)
#    yFitHM=(yFit.max()+yFit.min())/2
#    x0=np.argmin(yFit)
#    yFWHMList.append(x0+np.argmax(yFit[x0:]>yFitHM)-np.argmax(yFit[:x0]<yFitHM))
#camera.close()
#print(yFWHMList)
#temp=np.asarray(yFWHMList)
#print('FWHM',np.round(np.mean(temp),1),' +/- ', np.round(np.std(temp),1))
#
#
#plt.plot(y)
##plt.plot(ySmooth)
#plt.plot(xArr,yFit)
#plt.grid()
#plt.show()
#
#
#
##plt.imshow(img)
##plt.show()