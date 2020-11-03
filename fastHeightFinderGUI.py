from astropy.io import fits
import time
from DAQClass import DAQPin
import globalVariables as gv
from CameraClass import Camera
import sys
import numpy as np
import tkinter as tk
import matplotlib.pyplot as plt
import os
import scipy.special as sps
import scipy.interpolate as spi
import scipy.optimize as spo



class GUI:
    def __init__(self):
        self.settingsList=[] #list of tkinter field objects to save to a settings file
        self.x1Box=None #tkinter box object for image parameters
        self.x2Box=None #tkinter box object for image parameters
        self.y1Box=None #tkinter box object for image parameters
        self.y2Box=None #tkinter box object for image parameters
        self.voltOnResArr=None  #array to to hold voltage values to scan over near resonance
        self.voltPlotArr=None #array to be used in plotting
        self.camera=None  #to hold the camera opject
        self.galvoOut=None #galvo voltage control
        self.shutterOut=DAQPin(gv.shutterPin)  #shutter control for OP laser
        self.sigmaGuess=.025 #guess of value for sigma in volts
        self.gammaGuess=1e-3 #guess for value of gamma in volts
        self.window=tk.Tk()
        self.window.title("Fast Height Finder")
        self.window.geometry('800x600')

        lbl1=tk.Label(self.window, text='Profile Center (v)')
        lbl1.grid(column=0, row=0,sticky='E')

        self.v0Box=tk.Entry(self.window)
        self.v0Box.config(width=5)
        self.v0Box.grid(column=1, row=0, sticky='W')
        self.settingsList.append(self.v0Box)


        lbl32=tk.Label(self.window, text='bin size')
        lbl32.grid(column=2, row=0)

        self.binSizeBox=tk.Entry(self.window)
        self.binSizeBox.config(width=5)
        self.binSizeBox.grid(column=3, row=0, sticky='W')
        self.settingsList.append(self.binSizeBox)



        lbl3=tk.Label(self.window, text='FWHM (v)')
        lbl3.grid(column=0, row=1,sticky='E')

        self.fwhmBox=tk.Entry(self.window)
        self.fwhmBox.config(width=5)
        self.fwhmBox.grid(column=1, row=1, sticky='W')
        self.settingsList.append(self.fwhmBox)

        lbl3=tk.Label(self.window, text='num images off resonance')
        lbl3.grid(column=0, row=2)

        self.numImgOffResBox=tk.Entry(self.window)
        self.numImgOffResBox.config(width=5)
        self.numImgOffResBox.grid(column=1, row=2, sticky='W')
        self.settingsList.append(self.numImgOffResBox)

        lbl31=tk.Label(self.window, text='Exp time (ms)')
        lbl31.grid(column=2, row=1)

        self.expTimeBox=tk.Entry(self.window)
        self.expTimeBox.config(width=5)
        self.expTimeBox.grid(column=3, row=1, sticky='W')
        self.settingsList.append(self.expTimeBox)


        lbl3=tk.Label(self.window, text='num images on resonance')
        lbl3.grid(column=0, row=3)

        self.numImgOnResBox=tk.Entry(self.window)
        self.numImgOnResBox.config(width=5)
        self.numImgOnResBox.grid(column=1, row=3, sticky='W')
        self.settingsList.append(self.numImgOnResBox)



        lbl5=tk.Label(self.window, text='image x1')
        lbl5.grid(column=0, row=4)
        self.x1Box=tk.Entry(self.window)
        self.x1Box.config(width=5)
        self.x1Box.grid(column=1, row=4, sticky='W')
        self.settingsList.append(self.x1Box)

        lbl5=tk.Label(self.window, text='image x2')
        lbl5.grid(column=2, row=4)
        self.x2Box=tk.Entry(self.window)
        self.x2Box.config(width=5)
        self.x2Box.grid(column=3, row=4, sticky='W')
        self.settingsList.append(self.x2Box)

        lbl5=tk.Label(self.window, text='image y1')
        lbl5.grid(column=0, row=5)
        self.y1Box=tk.Entry(self.window)
        self.y1Box.config(width=5)
        self.y1Box.grid(column=1, row=5, sticky='W')
        self.settingsList.append(self.y1Box)

        lbl5=tk.Label(self.window, text='image y2')
        lbl5.grid(column=2, row=5)
        self.y2Box=tk.Entry(self.window)
        self.y2Box.config(width=5)
        self.y2Box.grid(column=3, row=5, sticky='W')
        self.settingsList.append(self.y2Box)



        lbl4=tk.Label(self.window, text='Camera')
        lbl4.grid(column=0, row=6)

        self.cameraVar=tk.StringVar(self.window)
        self.cameraVar.set("FAR")
        cameraChoice=["NEAR", "FAR"]
        CAMERA_MENU=tk.OptionMenu(self.window, self.cameraVar, *cameraChoice)
        CAMERA_MENU.grid(column=1, row=6, columnspan=2,sticky='W')
        self.settingsList.append(self.cameraVar)


        lbl4=tk.Label(self.window, text='Shutter')
        lbl4.grid(column=2, row=6)

        self.shutterVar=tk.StringVar(self.window)
        self.shutterVar.set("OPEN")
        shutterChoice=["CLOSED", "OPEN"]
        shutter_MENU=tk.OptionMenu(self.window, self.shutterVar, *shutterChoice,command=self.toggle_Shutter)
        shutter_MENU.grid(column=3, row=6, columnspan=2,sticky='W')
        self.settingsList.append(self.shutterVar)


        self.saveDataVar=tk.BooleanVar()
        saveDataCheckButton=tk.Checkbutton(self.window, text='save data', variable=self.saveDataVar)
        saveDataCheckButton.grid(column=1, row=7)
        self.settingsList.append(self.saveDataVar)


        self.ratioVar=tk.BooleanVar()
        ratioVarButton=tk.Checkbutton(self.window, text='ratio', variable=self.ratioVar)
        ratioVarButton.grid(column=1, row=8)
        self.settingsList.append(self.ratioVar)

        self.showPlotVar=tk.BooleanVar()
        showDataAnalysiButton=tk.Checkbutton(self.window, text='Show plot',
                                             variable=self.showPlotVar)
        showDataAnalysiButton.grid(column=1, row=9)
        self.settingsList.append(self.showPlotVar)






        lbl3=tk.Label(self.window, text='Run name')
        lbl3.grid(column=0, row=15)

        self.fileName=tk.Entry(self.window)
        self.fileName.config(width=60)
        self.fileName.grid(column=1, row=15, sticky='W', columnspan=40)
        self.settingsList.append(self.fileName)

        lbl3=tk.Label(self.window, text='Folder path')
        lbl3.grid(column=0, row=16)

        self.folderPath=tk.Entry(self.window)
        self.folderPath.config(width=60)
        self.folderPath.grid(column=1, row=16, sticky='W', columnspan=40)
        self.settingsList.append(self.folderPath)

        runButton=tk.Button(self.window, text='RUN', font=("Arial", 20), background="green", command=self.run)
        runButton.config(height=2, width=10)
        runButton.grid(column=0, row=17, columnspan=4, rowspan=4)

        coolCameraButton=tk.Button(self.window, text='cool camera', background="royal blue", command=self.cool_Camera)
        #coolCameraButton.config(height=2, width=10)
        coolCameraButton.grid(column=0, row=21, columnspan=4, rowspan=4)

        self.load_Settings()

        self.window.protocol("WM_DELETE_WINDOW", self.close_GUI)
        self.window.mainloop()
    def toggle_Shutter(self,x):
        if x=='CLOSED':
            self.close_Aperture()
        if x=='OPEN':
            self.open_Aperture()
    def cool_Camera(self):
        #cool the far field camera down
        if self.cameraVar.get()=="NEAR":
            print("YOU CAN'T COOL THE NEAR FIELD CAMERA")
            gv.warning_Sound()
        else:
            tempCamera=Camera(self.cameraVar.get(), 1000)  #the camera will cool down
            tempCamera.close()  #now close it. It will stay cool though
    def open_Aperture(self):
        #open the OP shutter
        self.shutterOut.write_Low()
        time.sleep(.01)

    def close_Aperture(self):
        #close the OP shutter
        self.shutterOut.write_High()
        time.sleep(.01)

    def close_GUI(self):
        self.shutterOut.close()
        self.save_Settings()
        self.window.destroy()
        sys.exit()
    def initialize_Camera(self):
        x1=int(self.x1Box.get())  #image region values
        y1=int(self.y1Box.get())  #image region values
        x2=int(self.x2Box.get())  #image region values
        y2=int(self.y2Box.get())  #image region values
        imageParams=[x1, x2, y1, y2]
        binSize=int(self.binSizeBox.get())  #binning value

        #initialize camera object. can be near or far field camera
        expTime=int(self.expTimeBox.get())  #exposure time, ms
        whichCam=self.cameraVar.get()  #which camera to use, 'FAR' or 'NEAR'
        self.camera=Camera(whichCam, expTime, imageParams=imageParams, bin=binSize)
    def close_Camera(self):
        self.camera.close()
        self.camera=None
    def initialize_Scan_And_Plot_Arrays(self):
        v0=float(self.v0Box.get()) #center value of transition from user
        df=float(self.fwhmBox.get()) #fwhm value from user
        offResFact=3 #go this many fwhm away from center for 'far' off resonance background
        numImagesOnRes=int(self.numImgOnResBox.get()) #number of images to take near the resonance
        numImagesOffRes=int(self.numImgOffResBox.get())
        self.voltOffResArr=np.linspace(v0-offResFact*df/2,v0-offResFact*df/2,num=numImagesOffRes) #voltage value to take images at 'far' off resonance
        self.voltOnResArr=np.linspace(v0-df,v0+df,num=numImagesOnRes) #array of voltages to take images of near the peak
        self.voltPlotArr=np.linspace(v0-offResFact*df,v0+offResFact*df,num=1000)#voltages to make plot with. This should
            #be dense and uniform so it looks good
    def run(self):
        self.save_Settings()
        self.initialize_Camera()
        if os.path.isdir(self.folderPath.get())==False:
            print('-----ERROR-----------')
            print('YOU HAVE ENTERED AN INVALID FOLDERPATH')
        if os.path.isfile(self.folderPath.get()+'\\'+self.fileName.get()+'.png')==True:
            print('--------------ERROR-------')
            print('THERE IS ALREADY A FILE WITH THAT NAME IN THAT FOLDER')
            gv.error_Sound()
            sys.exit()

        if self.ratioVar.get()==True:
            self.sweep_With_Shutter()
        else:
            self.sweep_Without_Shutter()
        self.camera.close()

    def sweep_With_Shutter(self):
        self.initialize_Scan_And_Plot_Arrays()
        self.initialize_Camera()

        voltList=[] #list to hold voltage values of corresponding images
        signalList1=[] #list for signal values for apeture open
        signalList2=[]  #list for signal values for apeture open




        gv.begin_Sound(noWait=True)#beep without waiting after the beep
        self.galvoOut=DAQPin(gv.galvoOutPin) #open the galvo control pin
        #take images 'far' off resonance. This scan is very close to together
        for volt in self.voltOffResArr:
            self.galvoOut.write(volt) #move galvo to new position
            voltList.append(volt)  #record the voltage value

            self.open_Aperture() #'turn on' the optical pumping
            img=self.camera.aquire_Image()
            signalList1.append(np.mean(img))

            self.close_Aperture() #'turn off' the optical pumping
            img=self.camera.aquire_Image()
            signalList2.append(np.mean(img))
        #now sweep around the peak near resonance
        for volt in self.voltOnResArr:
            self.galvoOut.write(volt) #move galvo to new position
            voltList.append(volt) #record the voltage value

            self.open_Aperture() #'turn on' the optical pumping
            img=self.camera.aquire_Image() #capture an image
            signalList1.append(np.mean(img)) #add the average of the pixels

            self.close_Aperture() #'turn off' the optical pumping
            img=self.camera.aquire_Image()#capture an image
            signalList2.append(np.mean(img))#add the average of the pixels


        self.galvoOut.close() #close and zero the pin
        self.open_Aperture() #open the shutter up again when done
        gv.finished_Sound(noWait=True) #beep without waiting after the beep

        #convert lists to arrays
        signalArr1=np.asarray(signalList1)
        signalArr2=np.asarray(signalList2)
        voltArr=np.asarray(voltList)

        #fit the data and get the optimal parameters and the error
        params1, perr1=self.fit_Data(voltArr, signalArr1)
        params2, perr2=self.fit_Data(voltArr, signalArr2)

        plt.close('all')
        plt.figure(figsize=(13,8))
        plt.plot(self.voltPlotArr, self.spectral_Profile(self.voltPlotArr, *params1), c='blue', label='fit, opened shutter')
        plt.plot(self.voltPlotArr, self.spectral_Profile(self.voltPlotArr, *params2), c='red', label='fit, closed shutter')
        plt.scatter(voltArr, signalArr1, label='data, opened shutter',c='blue')
        plt.scatter(voltArr, signalArr2, label='data, closed shutter',c='red',marker='x',s=100)

        v0=(params1[0]+params2[0])/2 #Get the center from the average of the two centers
        floor=(params1[1]+params2[1])/2 #get the floor from the average of the two floors
        #now use the v0 above for when the user runs again. This helps compensate for the laser drifting without the user
        #having to
        self.v0Box.delete(0,'end') #clear existing number
        self.v0Box.insert(0,str(np.round(v0,3))) #insert new number
        plt.axvline(x=v0, c='r', linestyle=':')
        plt.grid()
        plt.text(v0, floor, np.round(float(v0), 3)) #TODO: WHY IS THIS NOT WORKING ONLY HERE?

        plt.legend()
        ratio=np.round(params1[1]/params2[1],2)
        error=np.round(ratio*np.sqrt((perr1[1]/params1[1])**2+(perr2[1]/params2[1])**2),3)
        plt.suptitle('Ratio of open to close shutter height = '+str(ratio)+' +/- ' +str(error))
        plt.title("shutter open= "+str(np.round(params1[1],2))+' +/-'+str(np.round(perr1[1],1))+" . shutter closed= "
                  +str(np.round(params2[1],2))+' +/-'+str(np.round(perr2[1],1)))
        if self.saveDataVar.get()==True:
            plt.savefig(self.folderPath.get()+'\\'+self.fileName.get()+'Graph.png')
        if self.showPlotVar.get()==True:
            plt.show()
    def sweep_Without_Shutter(self):
        self.initialize_Scan_And_Plot_Arrays()
        self.initialize_Camera()

        voltList=[]
        signalList=[]


        gv.begin_Sound(noWait=True)#beep without waiting after the beep
        self.galvoOut=DAQPin(gv.galvoOutPin)
        #take images 'far' off resonance. This scan is very close to together
        for volt in self.voltOffResArr:
            self.galvoOut.write(volt) #move the galvo to a new voltage value
            voltList.append(volt) #record the voltage
            img=self.camera.aquire_Image() #capture image
            signalList.append(np.mean(img)) #add the average of the image's pixels


        for volt in self.voltOnResArr:
            self.galvoOut.write(volt) #move the galvo to a new voltage value
            voltList.append(volt) #record the voltage
            img=self.camera.aquire_Image() #capture image
            signalList.append(np.mean(img)) #add the average of the image's pixels
        self.galvoOut.close() #close and zero the galvo
        self.camera.close()
        self.open_Aperture() #open the apeture up when done
        gv.finished_Sound(noWait=True)#beep without waiting after the beep
        voltArr=np.asarray(voltList)
        signalArr=np.asarray(signalList)

        params,perr=self.fit_Data(voltArr,signalArr)

        plt.close('all')
        plt.plot(self.voltPlotArr, self.spectral_Profile(self.voltPlotArr,*params), c='orange', label='fit')
        plt.scatter(voltArr, signalArr, label='data')

        v0=params[0]
        floor=params[2]
        #now use the v0 above for when the user runs again. This helps compensate for the laser drifting without the user
        #having to
        self.v0Box.delete(0,'end') #clear existing number
        self.v0Box.insert(0,str(np.round(v0,3))) #insert new number
        plt.axvline(x=v0, c='r', linestyle=':')
        plt.text(v0, floor, np.round(float(v0), 3))
        plt.legend()
        plt.grid()
        plt.title('Height = '+str(np.round(params[1], 1))+'+/- '+str(np.round(perr[1], 1)))

        if self.saveDataVar.get()==True:
            plt.savefig(self.folderPath.get()+'\\'+self.fileName.get()+'Graph.png')
        if self.showPlotVar.get()==True:
            plt.show()
    def fit_Data(self,x,y):
        #fit the data to get the optimal parameters
        x0Guess=float(self.v0Box.get())
        aGuess=np.max(y)-np.min(y)
        bGuess=np.min(y)
        guess=[x0Guess,aGuess,bGuess,self.sigmaGuess,self.gammaGuess]
        print(guess)
        bounds=([-5, 0, 0, 0, 0], [5, 100000, 100000, .1, .1])
        if bGuess>bounds[1][2] or aGuess>bounds[1][1]:
            print('GUESS VALUES ARE LARGER THAN BOUNDS. SIGNAL STRENGTH IS VERY HIGH')
            gv.error_Sound()
            sys.exit()
        params,pcov=spo.curve_fit(self.spectral_Profile,x,y,p0=guess,bounds=bounds)
        perr=np.sqrt(np.diag(pcov))

        return params,perr

    def spectral_Profile(self,x,x0,a,b,sigma,gamma):
        v0=sps.voigt_profile(0,sigma,gamma)
        v=sps.voigt_profile(x-x0,sigma,gamma)
        return a*(v/v0)+b
    def save_Settings(self):
        file=open("fastHeightFinderGUI_Settings.txt", "w")
        for item in self.settingsList:
            file.write(str(item.get())+'\n')
        file.close()

    def load_Settings(self):
        try:
            file=open("fastHeightFinderGUI_Settings.txt", "r")
        except:
            print("NO SETTINGS FILE FOUND")
            return
        i=0
        for item in file:
            item=item.strip()
            if i>=len(self.settingsList):
                pass
            else:
                if isinstance(self.settingsList[i], tk.StringVar):
                    self.settingsList[i].set(item)
                elif isinstance(self.settingsList[i], tk.Entry):
                    self.settingsList[i].insert(0, item)
                elif isinstance(self.settingsList[i], tk.BooleanVar):
                    if item=='False' or item=='True':
                        self.settingsList[i].set(item)
            i+=1
        file.close()
gui=GUI()
