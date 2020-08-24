import os.path
import time
import globalVariables as gv
import matplotlib.pyplot as plt
from astropy.io import fits
#import time
#import datetime
#import matplotlib
import sys
from Sweeper import Sweeper
import numpy as np
#import matplotlib.pyplot as plt
import tkinter as tk
import MakeMHzScale
import DataAnalysis


class ExperimentGUI:
    def __init__(self):
        self.DAQDataArr=None#array to hold data read from DAQ board. each row is a sample of data like
                # [outputVoltage, Lithium reference chamber voltage]
        self.imageArrList=None #list to hold array of images. looks like [[imageN1,imageF1],[imageN2,imageF2],..]
            #where each image is a numpy array. So first image in each pair is for near, second for far
        self.settingsList=[] #to store the variables whos values will be saved so that the gui won't be opened
            #with all the boxes empty!
        self.window = tk.Tk()
        self.window.title("aquisition and analysis")
        self.window.geometry('800x600')

        # -----------------------------------------------------------
        # -------------data analysis section-------------------------
        # -----------------------------------------------------------

        # ---------------row 0-----------------------
        lbl1 = tk.Label(self.window, text='camera')
        lbl1.grid(column=0, row=0)

        self.cameraVarData = tk.StringVar(self.window)
        self.settingsList.append(self.cameraVarData)
        self.cameraVarData.set("NEAR")
        cameraChoice = ["NEAR", "FAR", 'BOTH']
        CAMERA_MENU = tk.OptionMenu(self.window, self.cameraVarData, *cameraChoice)
        CAMERA_MENU.grid(column=1, row=0,columnspan=2)

        # ----------------row 1-------------------
        lbl2 = tk.Label(self.window, text='Near camera', font=("Times", 15))
        lbl2.grid(column=1, row=1,columnspan=3)

        lbl3 = tk.Label(self.window, text='Far camera', font=("Times", 15))
        lbl3.grid(column=5, row=1,columnspan=3)
        # ----------------row 2------------------
        dimBoxWidth = 8

        lbl4 = tk.Label(self.window, text='x1')
        lbl4.grid(column=0, row=2, sticky='E')
        self.x1NearBox = tk.Entry(self.window)
        self.settingsList.append(self.x1NearBox)
        self.x1NearBox.config(width=dimBoxWidth)
        self.x1NearBox.grid(column=1, row=2, sticky='W')

        lbl5 = tk.Label(self.window, text='x2')
        lbl5.grid(column=2, row=2, sticky='E')
        self.x2NearBox = tk.Entry(self.window)
        self.settingsList.append(self.x2NearBox)
        self.x2NearBox.config(width=dimBoxWidth)
        self.x2NearBox.grid(column=3, row=2, sticky='W')


        lbl4 = tk.Label(self.window, text='x1')
        lbl4.grid(column=5, row=2, sticky='E')
        self.x1FarBox = tk.Entry(self.window)
        self.settingsList.append(self.x1FarBox)
        self.x1FarBox.config(width=dimBoxWidth)
        self.x1FarBox.grid(column=6, row=2, sticky='W')

        lbl4 = tk.Label(self.window, text='x2')
        lbl4.grid(column=7, row=2, sticky='E')
        self.x2FarBox = tk.Entry(self.window)
        self.settingsList.append(self.x2FarBox)
        self.x2FarBox.config(width=dimBoxWidth)
        self.x2FarBox.grid(column=8, row=2, sticky='W')



        # --------------row 3------------------------
        lbl5 = tk.Label(self.window, text='y1')
        lbl5.grid(column=0, row=3, sticky='E')
        self.y1NearBox = tk.Entry(self.window)
        self.settingsList.append(self.y1NearBox)
        self.y1NearBox.config(width=dimBoxWidth)
        self.y1NearBox.grid(column=1, row=3, sticky='W')

        lbl6 = tk.Label(self.window, text='y2')
        lbl6.grid(column=2, row=3, sticky='E')
        self.y2NearBox = tk.Entry(self.window)
        self.settingsList.append(self.y2NearBox)
        self.y2NearBox.config(width=dimBoxWidth)
        self.y2NearBox.grid(column=3, row=3, sticky='W')

        lbl7 = tk.Label(self.window, text='y1')
        lbl7.grid(column=5, row=3, sticky='E')
        self.y1FarBox = tk.Entry(self.window)
        self.settingsList.append(self.y1FarBox)
        self.y1FarBox.config(width=dimBoxWidth)
        self.y1FarBox.grid(column=6, row=3, sticky='W')

        lbl8 = tk.Label(self.window, text='y2')
        lbl8.grid(column=7, row=3, sticky='E')
        self.y2FarBox = tk.Entry(self.window)
        self.settingsList.append(self.y2FarBox)
        self.y2FarBox.config(width=dimBoxWidth)
        self.y2FarBox.grid(column=8, row=3, sticky='W')



        # -------------row 4-----------------
        lbl9 = tk.Label(self.window, text='exp time (ms)')
        lbl9.grid(column=0, row=4, sticky='W')

        self.expNearBox = tk.Entry(self.window)
        self.settingsList.append(self.expNearBox)
        self.expNearBox.config(width=dimBoxWidth)
        self.expNearBox.grid(column=1, row=4, sticky='W')

        lbl10 = tk.Label(self.window, text='exp time (ms)')
        lbl10.grid(column=5, row=4, sticky='W')

        self.expFarBox = tk.Entry(self.window)
        self.settingsList.append(self.expFarBox)
        self.expFarBox.config(width=dimBoxWidth)
        self.expFarBox.grid(column=6, row=4, sticky='W')
        #-----------------row 5-------------------------

        lbl9 = tk.Label(self.window, text='bin size')
        lbl9.grid(column=0, row=5, sticky='W')

        self.binSizeNearBox = tk.Entry(self.window)
        self.settingsList.append(self.binSizeNearBox)
        self.binSizeNearBox.config(width=dimBoxWidth)
        self.binSizeNearBox.grid(column=1, row=5, sticky='W')

        lbl10 = tk.Label(self.window, text='bin size')
        lbl10.grid(column=5, row=5, sticky='W')

        self.binSizeFarBox = tk.Entry(self.window)
        self.settingsList.append(self.binSizeFarBox)
        self.binSizeFarBox.config(width=dimBoxWidth)
        self.binSizeFarBox.grid(column=6, row=5, sticky='W')

        # ------------row 6--------------
        lbl10 = tk.Label(self.window, text='num exposure')
        lbl10.grid(column=0, row=6, sticky='W', columnspan=3)

        self.expNumBox = tk.Entry(self.window)
        self.settingsList.append(self.expNumBox)
        self.expNumBox.config(width=dimBoxWidth)
        self.expNumBox.grid(column=3, row=6, sticky='W')

        # -------------row 7-------------

        lbl11 = tk.Label(self.window, text='image scan range(v)')
        lbl11.grid(column=0, row=7, sticky='W', columnspan=2)

        self.startVoltBox = tk.Entry(self.window)
        self.settingsList.append(self.startVoltBox)
        self.startVoltBox.config(width=dimBoxWidth)
        self.startVoltBox.grid(column=3, row=7, sticky='W')
        lbl111 = tk.Label(self.window, text='to')
        lbl111.grid(column=4, row=7, sticky='W')
        self.stopVoltBox = tk.Entry(self.window)
        self.settingsList.append(self.stopVoltBox)
        self.stopVoltBox.config(width=dimBoxWidth)
        self.stopVoltBox.grid(column=5, row=7, sticky='W')

        
        # --------------row 8---------------
        lbl12 = tk.Label(self.window, text='folder path')
        lbl12.grid(column=0, row=8, sticky='W', columnspan=2)

        self.dataFolderPathBox = tk.Entry(self.window)
        self.settingsList.append(self.dataFolderPathBox)
        self.dataFolderPathBox.config(width=50)
        self.dataFolderPathBox.grid(column=2, row=8, sticky='W', columnspan=5)
        # ----------------row 9---------------------
        lbl13 = tk.Label(self.window, text='file name')
        lbl13.grid(column=0, row=9, sticky='W', columnspan=2)

        self.dataFileNameBox = tk.Entry(self.window)
        self.settingsList.append(self.dataFileNameBox)
        self.dataFileNameBox.config(width=20)
        self.dataFileNameBox.grid(column=2, row=9, sticky='W', columnspan=3)

        # -------------RUN BUTTON-------------------

        RUN_EXPERIMENT = tk.Button(self.window, text='RUN', font=("Arial", 20), background="green",command=self.aquire_Data)
        RUN_EXPERIMENT.config(height=2, width=10)
        RUN_EXPERIMENT.grid(column=0, row=11, columnspan=4, rowspan=4)

        # --------------------FLOW METER-------------------
        lbl14 = tk.Label(self.window, text='Desired Helium flow (SCCM)')
        lbl14.grid(column=10, row=2, sticky='W', columnspan=2)

        self.flowRateBox = tk.Entry(self.window)
        self.settingsList.append(self.flowRateBox)
        self.flowRateBox.config(width=dimBoxWidth)
        self.flowRateBox.grid(column=12, row=2, sticky='W')

        flowBtn = tk.Button(self.window, text='FLOW', command=self.make_Flow)
        flowBtn.grid(column=13, row=2, columnspan=2)

        self.flowRateDisplayVar = tk.StringVar()
        self.settingsList.append(self.flowRateDisplayVar)
        self.flowRateDisplayVar.set("0.0")

        lbl151 = tk.Label(self.window, text='Current flowrate (SCCM)')
        lbl151.grid(column=10, row=3, sticky='W', columnspan=2)
        lbl15 = tk.Label(self.window, textvariable=self.flowRateDisplayVar)
        lbl15.grid(column=12, row=3, sticky='W', columnspan=2)

        # ----------------lambdip check-----------

        btn1 = tk.Button(self.window, text='Li reference check', command=self.Li_Reference_Check)
        btn1.grid(column=0, row=15)

        #-----------------------------------------------------------
        #-------------data analysis section-------------------------
        #-----------------------------------------------------------
        lbl16=tk.Label(self.window, text="---------image analysis----------", font=('Arial', 20))
        lbl16.grid(column=0, row=20, sticky='W', columnspan=10)

        #----------------row 21------------
        self.cameraVarAnl=tk.StringVar(self.window)
        self.settingsList.append(self.cameraVarAnl)
        self.cameraVarAnl.set("NEAR")
        self.settingsList.append(self.cameraVarAnl)
        cameraChoiceAnl=["NEAR", "FAR"]
        CAMERA_MENU_Anl=tk.OptionMenu(self.window, self.cameraVarAnl, *cameraChoiceAnl)
        CAMERA_MENU_Anl.grid(column=0, row=21, columnspan=1)

        #---------------cheat method-----------------
        #cheatVar=tk.IntVar()
        #chk1=tk.Checkbutton(self.window, text="cheat", variable=cheatVar)
        #chk1.grid(column=2, row=21, columnspan=1)
#
        #lbl25=tk.Label(self.window, text='offset')
        #lbl25.grid(column=3, row=21, sticky='W', columnspan=1)
        #showFitVar=tk.IntVar()
        #chk2=tk.Checkbutton(self.window, text="dont show fit", variable=showFitVar)
        #chk2.grid(column=5, row=21, columnspan=1)
#
        #OFFSET_BOX=tk.Entry(self.window)
        #OFFSET_BOX.config(width=10)
        #OFFSET_BOX.grid(column=4, row=21, sticky='W', columnspan=1)

        #--------------row 22------------------------------


        lbl20=tk.Label(self.window, text='folder path')
        lbl20.grid(column=0, row=22, sticky='W', columnspan=2)

        self.anlFolderPathBox=tk.Entry(self.window)
        self.anlFolderPathBox.config(width=50)
        self.anlFolderPathBox.grid(column=2, row=22, sticky='W', columnspan=5)
        self.settingsList.append(self.anlFolderPathBox)

        #----------------row 23---------------------


        lbl21=tk.Label(self.window, text='file name')
        lbl21.grid(column=0, row=23, sticky='W', columnspan=2)

        self.anlFileNameBox=tk.Entry(self.window)
        self.anlFileNameBox.config(width=20)
        self.anlFileNameBox.grid(column=2, row=23, sticky='W', columnspan=2)


        #--------------row 24------------------------------

        lbl6=tk.Label(self.window, text='x1')
        lbl6.grid(column=0, row=24, sticky='E')
        self.X1Anl=tk.Entry(self.window)
        self.X1Anl.config(width=dimBoxWidth)
        self.X1Anl.grid(column=1, row=24, sticky='W')

        lbl7=tk.Label(self.window, text='x2')
        lbl7.grid(column=2, row=24, sticky='E')
        self.X2Anl=tk.Entry(self.window)
        self.X2Anl.config(width=dimBoxWidth)
        self.X2Anl.grid(column=3, row=24, sticky='W')

        #--------------------row 25------------
        lbl8=tk.Label(self.window, text='y1')
        lbl8.grid(column=0, row=25, sticky='E')
        self.Y1Anl=tk.Entry(self.window)
        self.Y1Anl.config(width=dimBoxWidth)
        self.Y1Anl.grid(column=1, row=25, sticky='W')

        lbl9=tk.Label(self.window, text='y2')
        lbl9.grid(column=2, row=25, sticky='E')
        self.Y2Anl=tk.Entry(self.window)
        self.Y2Anl.config(width=dimBoxWidth)
        self.Y2Anl.grid(column=3, row=25, sticky='W')

        #-------------row 26----------------
        btn4=tk.Button(self.window, text='analyze', font=("Arial", 10), background="orange", command=self.do_Data_Analysis)
        btn4.grid(column=0, row=26, columnspan=1)

        self.load_Settings()
        #it's a pain when the directory is empty because I can't remember which way the slashes go. If the directory
        #is empty go ahead and fill it with what I would almost certainly want to start with. Only works on windows
        if self.dataFolderPathBox.get()=='':
            self.dataFolderPathBox.insert(0,'C:\Data')
        if self.dataFileNameBox.get()=='':
            self.dataFileNameBox.insert(0,'sumFkngName')

        if self.anlFolderPathBox.get()=='':
            self.anlFolderPathBox.insert(0,'C:\Data')
        if self.anlFileNameBox.get()=='':
            self.anlFileNameBox.insert(0,'sumFkngName')

        self.window.protocol("WM_DELETE_WINDOW", self.close_GUI)
        self.update_Flow_Rate()
        self.window.mainloop()
    def update_Flow_Rate(self):
        #this function calls itself every 500 ms
        self.window.after(500, self.update_Flow_Rate)
    def close_GUI(self):
        self.save_Settings()
        #flowOut = DAQPin(gv.flowOutPin)  # want to make sure that flowrate is zero when program closes
        #flowOut.write(0.0)  # zero flowrate by writing zero voltage
        #flowOut.close()
        self.window.destroy()
        sys.exit()
    def catch_Errors(self):
        #this checks all the boxes that should be numbers to see if they contain a number. If this isn't done here
        #then a less helpful error will be thrown when trying to convert nonsense into a number
        tempList=[]
        tempList.append(self.x1NearBox.get())
        tempList.append(self.x2NearBox.get())
        tempList.append(self.y1NearBox.get())
        tempList.append(self.y2NearBox.get())
        tempList.append(self.x1FarBox.get())
        tempList.append(self.x2FarBox.get())
        tempList.append(self.y1FarBox.get())
        tempList.append(self.y2FarBox.get())
        tempList.append(self.expFarBox.get())
        tempList.append(self.expNearBox.get())
        tempList.append(self.expNumBox.get())
        tempList.append(self.binSizeFarBox.get())
        tempList.append(self.binSizeNearBox.get())
        tempList.append(self.startVoltBox.get())
        tempList.append(self.stopVoltBox.get())
        for item in tempList:
            if item=='':
                print('----------ERROR-------------')
                print('YOU HAVE NOTHING ENTERED FOR ONE OF THE CAMERA/IMAGE/SCAN PARAMETER BOXS')
                gv.error_Sound()
                sys.exit()
            try:
                int(item)
            except:
                print('----------ERROR-------------')
                print('YOU HAVE ENTERED SOMETHING OTHER THAN A NUMBER FOR ONE OF THE CAMERA/IMAGE/SCAN PARAMETER BOXS')
                print('ITS POSSIBLE THAT YOU ENTERED SOMETHING WITH A WEIRD FORMAT ALSO')
                gv.error_Sound()
                sys.exit()
        folderPath=self.dataFolderPathBox.get()
        fileName=self.dataFileNameBox.get()
        if folderPath=='':
            print('-----ERROR-----------')
            print('YOU HAVENT ENTERED ANYTHING FOR THE FOLDERPATH')
            gv.error_Sound()
            sys.exit()
        if fileName=='':
            print('-----ERROR-----------')
            print('YOU HAVENT ENTERED ANYTHING FOR THE FOLDERPATH')
            gv.error_Sound()
            sys.exit()

        if os.path.isdir(folderPath)==False:
            print('-----ERROR-----------')
            print('YOU HAVE ENTERED AN INVALID FOLDERPATH')
            gv.error_Sound()
            sys.exit()

        if os.path.isfile(folderPath+'\\'+fileName+'Near.fits') ==True:
            print('-----ERROR-----------')
            print('YOU ALREADY HAVE A FILE WITH THAT NAME IN THAT FOLDER')
            gv.error_Sound()
            sys.exit()
        if os.path.isfile(folderPath+'\\'+fileName+'Far.fits') ==True:
            print('-----ERROR-----------')
            print('YOU ALREADY HAVE A FILE WITH THAT NAME IN THAT FOLDER')
            gv.error_Sound()
            sys.exit()
    def do_Data_Analysis(self):
        pass
    def aquire_Data(self):
        self.catch_Errors()
        self.save_Settings()
        sweeper=Sweeper(self)
        sweeper.sweep()
        self.save_Fits_Files()
        self.save_Data()
        self.make_Info_File()
        #self.analyze_Data()
    def make_Info_File(self):
        return
    def save_Data(self):
        folderPath=self.dataFolderPathBox.get()
        fileName=self.dataFileNameBox.get()
        np.savetxt(folderPath+'\\'+fileName+'DAQData.csv',self.DAQDataArr,delimiter=',')
    def save_Fits_Files(self):
        folderPath=self.dataFolderPathBox.get()
        fileName=self.dataFileNameBox.get()
        nearImageList=[]
        farImageList=[]
        for item in self.imageArrList:
            #images need to be flipped to conform to FIT file conventions
            if item[0] is not None:
                nearImageList.append(np.rot90(np.transpose(item[0])))
            if item[1] is not None:
                farImageList.append(np.rot90(np.transpose(item[1])))
        if len(nearImageList)!=0:
            hdu=fits.PrimaryHDU(nearImageList)
            hdul=fits.HDUList([hdu])
            hdul.writeto(folderPath+'\\'+fileName+'Near.fits')
        if len(farImageList)!=0:
            hdu=fits.PrimaryHDU(farImageList)
            hdul=fits.HDUList([hdu])
            hdul.writeto(folderPath+'\\'+fileName+'Far.fits')
        print('images saved to fits files')
    def analyze_Data(self):
        MHzScaleArr=MakeMHzScale.make_MHz_Scale(self.DAQData)
        #PF=DataAnalysis.fitPixelData(MHzScaleArr)
        plt.plot(MHzScaleArr,self.DAQData[:,1])
        plt.show()


    def make_Flow(self):
        None
    def Li_Reference_Check(self):
        None
    def load_Settings(self):
        try:
            file=open("GUI_Settings.txt", "r")
        except:
            print("NO SETTINGS FILE FOUND")
            return
        i=0
        for item in file:
            item=item.strip()
            if i>=len(self.settingsList)-1:
                None
            else:
                if isinstance(self.settingsList[i],tk.StringVar):
                    self.settingsList[i].set(item)
                elif isinstance(self.settingsList[i],tk.Entry):
                    self.settingsList[i].insert(0, item)
            i+=1
        file.close()
    def save_Settings(self):
        file = open("GUI_Settings.txt", "w")
        for item in self.settingsList:
            file.write(item.get()+'\n')
        file.close()

experiment=ExperimentGUI()