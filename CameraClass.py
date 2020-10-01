import sys
import globalVariables as gv
import matplotlib.pyplot as plt
import pyfli as fli
import time
import progressbar
import numpy as np
from ximea import xiapi as xi
import threading



#TODO: DONT PRESENT A WARNING FOR LEAVING imageParams=None
#TODO: FIX PROGRESS BAR NOT WORKING

class Camera:
    def __init__(self,camName,expTime,imageParams=None,bin=None,temp=-25,binx=None,biny=None):
        #camName: Name of camera, or rather position. Far field vs near field
        #expTime: exposure time of camera in milliseconds
        #imageParam: A list of image paramters, [x1,x2,y1,y2], where x1 is start of image region along x, x2 is end of
                #image region along x. y1 is start of image region along y, y2 is end of image region along y.
                #if none then use the entire sensor
        #bin: Bin size in both directions, ie square binning. For Ximea this is done to the numpy array output. For the
                #FLI camera it is done with hardware on the camera. For a CMOS camera hardware binning doesn't improve
                #performance like it does on a CCD. If the value is None then it set to 1 in the catch errors section
        #binx: binning in the x direction
        #biny: binning in the y direction

        self.camName=camName
        self.expTime=expTime
        if imageParams==None:
            imageParams=[0,9999999,0,999999] #set to full range
        self.imageParams=imageParams
        self.bin=bin
        self.binx=binx
        self.biny=biny
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
        if self.imageParams is None:
            #if no value is provided use the full feild of view
            x2=self.cameraObject.get_width_maximum()
            y2=self.cameraObject.get_height_maximum()
            self.imageParams=[0,x2,0,y2]

        self._catch_And_Fix_Errors()
        self.cameraObject.set_imgdataformat('XI_MONO16')  #return data with the greatest depth possible. Otherwise data
        #is forced to fit into an 8 bit array
        self.cameraObject.set_exposure(int(self.expTime*1E3)) #exposure time is in microseconds for ximea camera
        #Turn off every feature that I can. Who knows whats happening if these are on? I'm not even confident in this
        #camera with these off
        self.cameraObject.set_gain(0)
        self.cameraObject.disable_bpc()
        self.cameraObject.disable_ffc()
        self.cameraObject.disable_aeag()
        self.cameraObject.disable_LUTEnable()
        #----disable obnoxious image clamping thing----
        self.cameraObject.set_param('device_unit_register_selector', 0x0222)
        self.cameraObject.set_param('device_unit_register_value', 0xF0)

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
        if self.camName!='FAR' and self.camName!='NEAR':
            print('---------------ERROR----------------')
            print('NO VALID CAMERA NAME PROVIDED. VALID NAMES ARE \'NEAR\' AND \'FAR\'')
            gv.error_Sound()
            sys.exit()


        #binning results cannot conflict
        binError=False
        if self.bin is None and self.binx is None and self.biny is None:
            binError=True
        if self.bin is not None and (self.binx is not None or self.biny is not None):
            binError=True
        if self.bin is None and (self.binx is None or self.biny is None):
            binError=True
        if binError==True:
            print('------------ERROR-------------')
            print('YOU HAVE SET A DISALLOWED BIN CONFIGURATION')
            print('ALLOWED CONFIGURATIONS: bin=value, binx=None,biny=None OR bin=None, binx=value1,biny=value2')
            print('BY DEFAULT ALL ARE NONE AND THEN bin IS THEN SET TO 1')
            sys.exit()
        if self.bin is None and self.binx is None and self.biny is None:
            self.bin=1
            self.binx=1
            self.biny=1
        if self.bin is not None:
            self.binx=self.bin
            self.biny=self.bin

        x1, x2, y1, y2=self.imageParams
        if self.camName=='FAR':
            if self.binx<0 or self.binx>16 or self.biny<0 or self.biny>16:
                print('----------ERROR-------------')
                print('VALID RANGE OF BINNING FOR FLI CAMERA IS 1 TO 16')
                gv.error_Sound()
                sys.exit()
            if self.expTime<50:
                print('----------ERROR-------------')
                print('FLI CAMERA EXPOSURE MUST BE AT LEAST 50 ms')
                gv.error_Sound()
                sys.exit()


        if self.camName=='NEAR':
            xMax=self.cameraObject.get_width_maximum()
            yMax=self.cameraObject.get_height_maximum()
        if self.camName=='FAR':
            xMax=1024
            yMax=1024

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
            x1New=self.binx*8*(x1//(self.binx*8))
            x2New=self.binx*8*(x2//(self.binx*8))
            y1New=self.biny*8*(y1//(self.biny*8))
            y2New=self.biny*8*(y2//(self.biny*8))
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
            x1New=self.binx*(x1//self.binx)
            x2New=self.binx*(x2//self.binx)
            y1New=self.biny*(y1//self.biny)
            y2New=self.biny*(y2//self.biny)
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

        try:
            self.cameraIndex=fli.FLIOpen('flipro0', 'usb', 'camera')
        except:
            raise Exception('NO FLI CAMERA CONNECTED, OR INCORRECT FLI CAMERA CONNECTED. READ COMMENTS BELOW IN CODE FOR MORE INFO')
            #For some reason when trying to get a list of connected FLI cameras and their names if there are none
            #connected then a silent error occurs and python exits. This does not occur if I do FLIOpen first, but it
            #depends on me knowing the name of the camera which I think is always 'flipro0' for the first one. This error
            #will be thrown if that isn't true as well, such as using a different FLI camera, or some other change
            #in ports or something. Wish it didn't have to be like this but this library is a home made wrapper
        if self.imageParams is None:
            #if no value is provided use the full feild of view
            self.imageParams=[0,1024,0,1024]
        self._catch_And_Fix_Errors()
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
        fli.setVBin(self.cameraIndex,self.binx)
        fli.setHBin(self.cameraIndex,self.biny)
        fli.controlBackgroundFlush(self.cameraIndex,'start')

        #This part is a little tricky. The image is rotated 90 degrees from the camera mounting bracket so I transform
        #the coordinates into the rotated system.
        x1,x2,y1,y2=self.imageParams #coordinates in upright frame

        #Rotate the coordinates. Because the frame is specified from the bottom left like in a fits file instead of
        #the top left like the FLI camera wants I need to adjust for that in regards to the y dimension. Numpy has the
        #same convention as the FLI camera. I have to do this in two steps to prevent variables from being overwritten
        #because obviously I can't do something like x1=y1, x2=y2, y1=x2,y2=x1. The best way to see this logic is to
        #draw a box with coordinates and see how it changes.
        x1Rot=y1
        x2Rot=y2
        y1Rot=1024-x2
        y2Rot=1024-x1
        binx=self.biny #binning dimension flips
        biny=self.binx

        x1=x1Rot
        x2=x2Rot
        y1=y1Rot
        y2=y2Rot
        y1Temp=1024-y2
        y2Temp=1024-y1
        y1=y1Temp
        y2=y2Temp

        #now adjust to conform to the input dimensions for the FLI camera
        x2=x1+(x2-x1)//binx
        y2=y1+(y2-y1)//biny

        #some of this could be streamlined, like removing the redundant 1024 that gets subtracted anyways, but I want
        #to preserve the flow of the logic
        fli.setImageArea(self.cameraIndex, x1, y1, x2, y2) #set the image area
    def aquire_Image(self):
        if self.camName=='FAR':
            img=self._aquire_Image_FLI()
            img=np.rot90(img) #image is sideways relative to mounting screw, needs to be rotated
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
        img=self._bin_Image(img)
        return img
    def _aquire_Image_FLI(self):
        fli.exposeFrame(self.cameraIndex)
        while fli.getExposureStatus(self.cameraIndex)>0: #when there is 0 milliseconds remaining the exposure is done
            time.sleep(.01) #I don't think I want to query the camera for exposure status millions of time a second
        img=fli.grabFrame(self.cameraIndex)
        return img

    def _bin_Image(self,image):
        #ximea images have to be manually binned unfortunately...
        m, n=image.shape
        return image.reshape(m//self.biny, self.biny, n//self.binx, self.binx).sum(3).sum(1)
    def close(self):
        if self.camName=='NEAR':
            self.cameraObject.stop_acquisition()
            self.cameraObject.close_device()
