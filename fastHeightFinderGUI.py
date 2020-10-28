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
        self.settingsList=[]
        self.x1Box=None
        self.x2Box=None
        self.y1Box=None
        self.y2=None
        self.voltOnResArr=None  #array to to hold voltage values to scan over near resonance
        self.voltOffRes=None #voltage value far off resonance
        self.voltPlotArr=None #array to be used in plotting
        self.camera=None  #to hold the camera opject
        self.galvoOut=None
        self.shutterOut=None  #shutter control
        self.sigmaGuess=.025 #guess of value for sigma in volts
        self.gammaGuess=1e-3 #guess for value of gamma in volts
        self.window=tk.Tk()
        self.window.title("Simple Scan")
        self.window.geometry('800x600')

        lbl1=tk.Label(self.window, text='Profile Center (V)')
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



        lbl3=tk.Label(self.window, text='FWHM')
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


        self.saveDataVar=tk.BooleanVar()
        saveDataCheckButton=tk.Checkbutton(self.window, text='save data', variable=self.saveDataVar)
        saveDataCheckButton.grid(column=1, row=7)
        self.settingsList.append(self.saveDataVar)


        self.shutterVar=tk.BooleanVar()
        shutterVarButton=tk.Checkbutton(self.window, text='shutter', variable=self.shutterVar)
        shutterVarButton.grid(column=1, row=8)
        self.settingsList.append(self.shutterVar)

        self.showPlotVar=tk.BooleanVar()
        showDataAnalysiButton=tk.Checkbutton(self.window, text='Show plot',
                                             variable=self.showPlotVar)
        showDataAnalysiButton.grid(column=1, row=9)
        self.settingsList.append(self.showPlotVar)






        lbl3=tk.Label(self.window, text='Run name')
        lbl3.grid(column=0, row=15)

        self.fileName=tk.Entry(self.window)
        self.fileName.config(width=20)
        self.fileName.grid(column=1, row=15, sticky='W', columnspan=20)
        self.settingsList.append(self.fileName)

        lbl3=tk.Label(self.window, text='Folder path')
        lbl3.grid(column=0, row=16)

        self.folderPath=tk.Entry(self.window)
        self.folderPath.config(width=30)
        self.folderPath.grid(column=1, row=16, sticky='W', columnspan=30)
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

    def cool_Camera(self):
        if self.cameraVar.get()=="NEAR":
            print("YOU CAN'T COOL THE NEAR FIELD CAMERA")
            gv.warning_Sound()
        else:
            tempCamera=Camera(self.cameraVar.get(), 1000)  #the camera will cool down
            tempCamera.close()  #now close it. It will stay cool though

    def open_Aperture(self):
        self.shutterOut.write_Low()
        time.sleep(.05)

    def close_Aperture(self):
        self.shutterOut.write_High()
        time.sleep(.05)

    def close_GUI(self):
        self.save_Settings()
        self.window.destroy()
        sys.exit()

    def run(self):
        self.save_Settings()
        self.galvoOut=DAQPin(gv.galvoOutPin)
        self.shutterOut=DAQPin(gv.shutterPin)
        v0=float(self.v0Box.get()) #center value of transition from user
        df=float(self.fwhmBox.get()) #fwhm value from user
        offResFact=3 #go this many fwhm away from center
        numOnRes=int(self.numImgOnResBox.get()) #number of images to take near the resonance
        x1=int(self.x1Box.get())
        y1=int(self.y1Box.get())
        x2=int(self.x2Box.get())
        y2=int(self.y2Box.get())
        imageParams=[x1,x2,y1,y2]
        binSize=int(self.binSizeBox.get())

        self.camera=Camera(self.cameraVar.get(), float(self.expTimeBox.get()), imageParams=imageParams, bin=binSize)

        self.voltOffRes=v0-offResFact*df #voltage value to take images at 'far' off resonance
        self.voltOnResArr=np.linspace(v0-df,v0+df,num=numOnRes) #array of voltages to take images of near the peak

        self.voltPlotArr=np.linspace(v0-offResFact*df,v0+offResFact*df,num=1000)#voltages to make plot with. This should
            #be dense and uniform so it looks good
        if os.path.isdir(self.folderPath.get())==False:
            print('-----ERROR-----------')
            print('YOU HAVE ENTERED AN INVALID FOLDERPATH')
        if os.path.isfile(self.folderPath.get()+'\\'+self.fileName.get()+'.png')==True:
            print('--------------ERROR-------')
            print('THERE IS ALREADY A FILE WITH THAT NAME IN THAT FOLDER')
            gv.error_Sound()
            sys.exit()
        plt.close('all')
        if self.shutterVar.get()==True:
            self.sweep_With_Shutter()
        else:
            self.sweep_Without_Shutter()
    def sweep_With_Shutter(self):

        voltList=[] #list to hold voltage values of corresponding images
        signalList1=[] #list for signal values for apeture open
        signalList2=[]  #list for signal values for apeture open



        #take images far off resonance. I don't waste time sweeping the laser here, just pile the values up at the same
        #voltage.
        self.galvoOut.write(self.voltOffRes)
        numImagesFarOff=int(self.numImgOffResBox.get())
        for i in range(numImagesFarOff):
            self.open_Aperture()
            img=self.camera.aquire_Image()
            signalList1.append(np.mean(img))

            self.close_Aperture()
            img=self.camera.aquire_Image()
            signalList2.append(np.mean(img))

            voltList.append(self.voltOffRes)


        #now sweep around the peak
        for volt in self.voltOnResArr:
            self.galvoOut.write(volt)
            img=self.camera.aquire_Image()
            self.open_Aperture()
            signalList1.append(np.mean(img))

            self.close_Aperture()
            img=self.camera.aquire_Image()
            signalList2.append(np.mean(img))

            voltList.append(volt)
        gv.finished_Sound()
        self.galvoOut.close()
        self.camera.close()

        signalArr1=np.asarray(signalList1)
        signalArr2=np.asarray(signalList2)
        voltArr=np.asarray(voltList)

        params1, perr1=self.fit_Data(voltArr, signalArr1)
        params2, perr2=self.fit_Data(voltArr, signalArr2)

        plt.figure(figsize=(13,8))
        plt.plot(self.voltPlotArr, self.spectral_Profile(self.voltPlotArr, *params1), c='blue', label='fit, opened shutter')
        plt.plot(self.voltPlotArr, self.spectral_Profile(self.voltPlotArr, *params2), c='red', label='fit, closed shutter')
        plt.scatter(voltArr, signalArr1, label='data, opened shutter',c='blue')
        plt.scatter(voltArr, signalArr2, label='data, closed shutter',c='red',marker='x',s=100)

        x0=(params1[0]+params2[0])/2
        floor=(params1[1]+params2[1])/2
        plt.axvline(x=x0, c='r', linestyle=':')
        plt.text(x0, floor, np.round(float(x0), 2))  #,transform=ax.transAxes)

        plt.legend()
        ratio=np.round(params1[1]/params2[1],2)
        plt.title('Ratio of open to close shutter height = '+str(ratio))
        if self.saveDataVar.get()==True:
            plt.savefig(self.folderPath.get()+'\\'+self.fileName.get()+'Graph.png')
        if self.showPlotVar.get()==True:
            plt.show()
    def sweep_Without_Shutter(self):
        voltList=[]
        signalList=[]



        #take images far off resonance. I don't waste time sweeping the laser here, just pile the values up at the same
        #voltage.
        self.galvoOut.write(self.voltOffRes)
        numImagesFarOff=int(self.numImgOffResBox.get()) #number of images far off resonance
        for i in range(numImagesFarOff):
            img=self.camera.aquire_Image()
            signalList.append(np.mean(img))
            voltList.append(self.voltOffRes)

        for volt in self.voltOnResArr:
            self.galvoOut.write(volt)
            voltList.append(volt)
            img=self.camera.aquire_Image()
            signalList.append(np.mean(img))
        gv.finished_Sound()
        self.galvoOut.close()
        self.camera.close()
        voltArr=np.asarray(voltList)
        signalArr=np.asarray(signalList)

        params,perr=self.fit_Data(voltArr,signalArr)

        plt.plot(self.voltPlotArr, self.spectral_Profile(self.voltPlotArr,*params), c='orange', label='fit')
        plt.scatter(voltArr, signalArr, label='data')

        x0=params[0]
        floor=params[2]
        plt.axvline(x=x0, c='r', linestyle=':')
        plt.text(x0, floor, np.round(float(x0), 2))  #,transform=ax.transAxes)
        plt.legend()
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
        bounds=([-5, 0, 0, 0, 0], [5, 10000, 10000, .1, .1])

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
