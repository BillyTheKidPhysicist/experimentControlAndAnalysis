import time
from DAQClass import DAQPin
import globalVariables as gv
from CameraClass import Camera
import sys
import numpy as np
import tkinter as tk
import matplotlib.pyplot as plt
import os


class GUI:
    def __init__(self):
        self.settingsList=[]
        self.x1Box=None
        self.x2Box=None
        self.y1Box=None
        self.y2=None
        self.voltArr=None  #array to to hold voltage value sto scan over
        self.camera=None #to hold the camera opject
        self.galvoOut=None
        self.shutterOut=None #shutter control
        self.window=tk.Tk()
        self.window.title("Simple Scan")
        self.window.geometry('800x600')

        lbl1=tk.Label(self.window, text='Scan Range (V)')
        lbl1.grid(column=0, row=0)

        self.voltStartBox=tk.Entry(self.window)
        self.voltStartBox.config(width=5)
        self.voltStartBox.grid(column=1, row=0, sticky='W')
        self.settingsList.append(self.voltStartBox)

        lbl2=tk.Label(self.window, text='to')
        lbl2.grid(column=2, row=0)

        self.voltStopBox=tk.Entry(self.window)
        self.voltStopBox.config(width=5)
        self.voltStopBox.grid(column=3, row=0, sticky='W')
        self.settingsList.append(self.voltStopBox)

        lbl3=tk.Label(self.window, text='Num images')
        lbl3.grid(column=0, row=1)

        self.numImagesBox=tk.Entry(self.window)
        self.numImagesBox.config(width=5)
        self.numImagesBox.grid(column=1, row=1, sticky='W')
        self.settingsList.append(self.numImagesBox)

        lbl31=tk.Label(self.window, text='Exp time (ms)')
        lbl31.grid(column=2, row=1)

        self.expTimeBox=tk.Entry(self.window)
        self.expTimeBox.config(width=5)
        self.expTimeBox.grid(column=3, row=1, sticky='W')
        self.settingsList.append(self.expTimeBox)

        lbl32=tk.Label(self.window, text='bin size')
        lbl32.grid(column=4, row=1)

        self.binSizeBox=tk.Entry(self.window)
        self.binSizeBox.config(width=5)
        self.binSizeBox.grid(column=5, row=1, sticky='W')
        self.settingsList.append(self.binSizeBox)

        lbl4=tk.Label(self.window, text='Camera')
        lbl4.grid(column=0, row=4)

        self.cameraVar=tk.StringVar(self.window)
        self.cameraVar.set("FAR")
        cameraChoice=["NEAR", "FAR"]
        CAMERA_MENU=tk.OptionMenu(self.window, self.cameraVar, *cameraChoice)
        CAMERA_MENU.grid(column=1, row=4, columnspan=2)
        self.settingsList.append(self.cameraVar)

        lbl5=tk.Label(self.window, text='image x1')
        lbl5.grid(column=0, row=5)
        self.x1Box=tk.Entry(self.window)
        self.x1Box.config(width=5)
        self.x1Box.grid(column=1, row=5, sticky='W')
        self.settingsList.append(self.x1Box)

        lbl5=tk.Label(self.window, text='image x2')
        lbl5.grid(column=2, row=5)
        self.x2Box=tk.Entry(self.window)
        self.x2Box.config(width=5)
        self.x2Box.grid(column=3, row=5, sticky='W')
        self.settingsList.append(self.x2Box)

        lbl5=tk.Label(self.window, text='image y1')
        lbl5.grid(column=0, row=6)
        self.y1Box=tk.Entry(self.window)
        self.y1Box.config(width=5)
        self.y1Box.grid(column=1, row=6, sticky='W')
        self.settingsList.append(self.y1Box)

        lbl5=tk.Label(self.window, text='image y2')
        lbl5.grid(column=2, row=6)
        self.y2Box=tk.Entry(self.window)
        self.y2Box.config(width=5)
        self.y2Box.grid(column=3, row=6, sticky='W')
        self.settingsList.append(self.y2Box)

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
        self.fileName.grid(column=1, row=15, sticky='W',columnspan=20)
        self.settingsList.append(self.fileName)

        lbl3=tk.Label(self.window, text='Folder path')
        lbl3.grid(column=0, row=16)

        self.folderPath=tk.Entry(self.window)
        self.folderPath.config(width=30)
        self.folderPath.grid(column=1, row=16, sticky='W',columnspan=30)
        self.settingsList.append(self.folderPath)

        runButton=tk.Button(self.window, text='RUN', font=("Arial", 20), background="green", command=self.run)
        runButton.config(height=2, width=10)
        runButton.grid(column=0, row=17, columnspan=4, rowspan=4)


        coolCameraButton=tk.Button(self.window, text='cool camera', background="royal blue",command=self.cool_Camera)
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
            tempCamera=Camera(self.cameraVar.get(),1000) #the camera will cool down
            tempCamera.close() #now close it. It will stay cool though
    def open_Aperture(self):
        self.shutterOut.write_Low()
        time.sleep(.001)
    def close_Aperture(self):
        self.shutterOut.write_High()
        time.sleep(.001)
    def close_GUI(self):
        self.save_Settings()
        self.window.destroy()
        sys.exit()

    def run(self):
        self.save_Settings()
        self.galvoOut=DAQPin(gv.galvoOutPin)
        self.shutterOut=DAQPin(gv.shutterPin)
        self.galvoOut.write(float(self.voltStartBox.get()))
        x1=float(self.x1Box.get())
        y1=float(self.y1Box.get())
        x2=float(self.x2Box.get())
        y2=float(self.y2Box.get())
        imageParams=[x1,x2,y1,y2]
        self.camera=Camera(self.cameraVar.get(),float(self.expTimeBox.get()),imageParams=imageParams)
        
        
        self.voltArr=np.linspace(int(self.voltStartBox.get()), int(self.voltStopBox.get()),
                                 num=int(self.numImagesBox.get()))
        if os.path.isdir(self.folderPath.get())==False:
            print('-----ERROR-----------')
            print('YOU HAVE ENTERED AN INVALID FOLDERPATH')
        if os.path.isfile(self.folderPath.get()+'\\'+self.fileName.get()+'.png')==True:
            print('--------------ERROR-------')
            print('THERE IS ALREADY A FILE WITH THAT NAME IN THAT FOLDER')
            gv.error_Sound()
            sys.exit()



        if self.shutterVar.get()==True:
            self.sweep_With_shutter()
        else:
            self.sweep_Without_shutter()

    def take_Dark_Image_Average(self, num=3):
        image=self.camera.aquire_Image()
        for i in range(num-1):
            image+=self.camera.aquire_Image()
        image=image/num  #average the three images
        return image

    def sweep_With_shutter(self):


        self.open_Aperture()
        darkImage1=self.take_Dark_Image_Average() #dark image shutter open
        self.close_Aperture()
        darkImage2=self.take_Dark_Image_Average() #dark image shutter closed

        y1List=[] #shutter open list of image sums
        y2List=[] #shutter closed list of image sums

        for volt in self.voltArr:
            self.galvoOut.write(volt)
            #take with shutter open
            self.open_Aperture()
            #gv.begin_Sound()
            image1=self.camera.aquire_Image()
            y1List.append(np.sum(image1-darkImage1))

            #take with shutter closed
            self.close_Aperture()
            #gv.finished_Sound()
            image2=self.camera.aquire_Image()
            y2List.append(np.sum(image2-darkImage2))
        self.galvoOut.close()
        self.shutterOut.close()
        self.camera.close()

        y1=np.asarray(y1List) #shutter open
        y2=np.asarray(y2List) #shutter closed
        ratio=np.trapz(y1)/np.trapz(y2)
        plt.suptitle('Signal with shutter closed and open')
        plt.title('ratio of integral of open to close = '+str(np.round(ratio, 3)))
        plt.plot(self.voltArr, y1, label='shutter open')
        plt.plot(self.voltArr, y2, label='shutter closed')
        plt.xlabel('Volts')
        plt.ylabel('Pixel counts')
        plt.legend()
        plt.grid()

        if self.saveDataVar.get()==True:
            plt.savefig(self.folderPath.get()+'\\'+self.fileName.get())
        if self.showPlotVar.get()==True:
            plt.show()
    def sweep_Without_shutter(self):
        darkImage=self.take_Dark_Image_Average()

        imageSumList=[]
        for volt in self.voltArr:
            self.galvoOut.write(volt)
            image=self.camera.aquire_Image()
            image=image-darkImage
            imageSumList.append(np.sum(image))
        imageSumArr=np.asarray(imageSumList)

        self.galvoOut.close()
        self.shutterOut.close()
        self.camera.close()

        plt.plot(self.voltArr, imageSumArr)
        plt.xlabel('Volts')
        plt.ylabel('Pixel counts')
        plt.legend()
        plt.grid()
        if self.saveDataVar.get()==True:
            plt.savefig(self.folderPath.get()+'\\'+self.fileName.get())
        if self.showPlotVar.get()==True:
            plt.show()
    def save_Settings(self):
        file=open("simpleScaneGUI_Settings.txt", "w")
        for item in self.settingsList:
            file.write(str(item.get())+'\n')
        file.close()

    def load_Settings(self):
        try:
            file=open("simpleScaneGUI_Settings.txt", "r")
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
#if os.path.isdir(folderPath) == False:
#    print('-----ERROR-----------')
#    print('YOU HAVE ENTERED AN INVALID FOLDERPATH')
#    gv.error_Sound()
#    sys.exit()