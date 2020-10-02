from datetime import datetime
import os.path
import time
from DAQClass import DAQPin
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
from Analyzer import Analyzer
import scipy.interpolate as spi


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

        self.expTimeNearBox = tk.Entry(self.window)
        self.settingsList.append(self.expTimeNearBox)
        self.expTimeNearBox.config(width=dimBoxWidth)
        self.expTimeNearBox.grid(column=1, row=4, sticky='W')

        lbl10 = tk.Label(self.window, text='exp time (ms)')
        lbl10.grid(column=5, row=4, sticky='W')

        self.expTimeFarBox = tk.Entry(self.window)
        self.settingsList.append(self.expTimeFarBox)
        self.expTimeFarBox.config(width=dimBoxWidth)
        self.expTimeFarBox.grid(column=6, row=4, sticky='W')
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
        cameraChoiceAnl=["NEAR", "FAR"]
        CAMERA_MENU_Anl=tk.OptionMenu(self.window, self.cameraVarAnl, *cameraChoiceAnl)
        CAMERA_MENU_Anl.grid(column=0, row=21, columnspan=1)

        lbl17=tk.Label(self.window, text="I/I_sat")
        lbl17.grid(column=1, row=21, sticky='E')

        self.saturationConstantBox=tk.Entry(self.window)
        self.saturationConstantBox.config(width=dimBoxWidth)
        self.saturationConstantBox.grid(column=2, row=21, sticky='W')
        self.settingsList.append(self.saturationConstantBox)

        #---------------cheat method-----------------
        #cheatVar=tk.IntVar()
        #cheatVar=tk.IntVar()file
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
        self.settingsList.append(self.anlFileNameBox)


        #--------------row 24------------------------------

        lbl6=tk.Label(self.window, text='x1')
        lbl6.grid(column=0, row=24, sticky='E')
        self.x1BoxAnl=tk.Entry(self.window)
        self.x1BoxAnl.config(width=dimBoxWidth)
        self.x1BoxAnl.grid(column=1, row=24, sticky='W')
        self.settingsList.append(self.x1BoxAnl)

        lbl7=tk.Label(self.window, text='x2')
        lbl7.grid(column=2, row=24, sticky='E')
        self.x2BoxAnl=tk.Entry(self.window)
        self.x2BoxAnl.config(width=dimBoxWidth)
        self.x2BoxAnl.grid(column=3, row=24, sticky='W')
        self.settingsList.append(self.x2BoxAnl)

        #--------------------row 25------------
        lbl8=tk.Label(self.window, text='y1')
        lbl8.grid(column=0, row=25, sticky='E')
        self.y1BoxAnl=tk.Entry(self.window)
        self.y1BoxAnl.config(width=dimBoxWidth)
        self.y1BoxAnl.grid(column=1, row=25, sticky='W')
        self.settingsList.append(self.y1BoxAnl)

        lbl9=tk.Label(self.window, text='y2')
        lbl9.grid(column=2, row=25, sticky='E')
        self.y2BoxAnl=tk.Entry(self.window)
        self.y2BoxAnl.config(width=dimBoxWidth)
        self.y2BoxAnl.grid(column=3, row=25, sticky='W')
        self.settingsList.append(self.y2BoxAnl)

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
        self.flowRateDisplayVar.set(str(self.get_Flow()))

    def get_Flow(self):
        flowIn=DAQPin(gv.flowInPin)
        b=4.29
        m=1.4236
        volt=flowIn.read()
        flowIn.close()
        return b+2000.0*(volt/5.0)*m
    def close_GUI(self):
        self.save_Settings()
        #flowOut = DAQPin(gv.flowOutPin)  # want to make sure that flowrate is zero when program closes
        #flowOut.write(0.0)  # zero flowrate by writing zero voltage
        #flowOut.close()
        self.window.destroy()
        sys.exit()
    def catch_Errors_Aquisition(self):
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
        tempList.append(self.expTimeFarBox.get())
        tempList.append(self.expTimeNearBox.get())
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
                float(item)
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
    def catch_Errors_Analysis(self):
        fileName=self.anlFileNameBox.get()
        folderPath=self.anlFolderPathBox.get()
        try:
            x1=int(self.x1BoxAnl.get())
            x2=int(self.x2BoxAnl.get())
            y1=int(self.y1BoxAnl.get())
            y2=int(self.y2BoxAnl.get())
        except:
            print('---------ERROR-----------')
            print('YOU MUST ENTER A NON NEGATIVE NUMBER IN THE IMAGE DIMENSION BOX.')
            print('YOU MAY HAVE LEFT THEM BLANK OR ENTER A LETTER, A SPACE, OR SOME OTHER NON NUMBER CHARACTER')
            gv.error_Sound()
            sys.exit()
        if x2<x1 or y2<y1:
            print('---------ERROR-------------')
            print('X2 AND Y2 MUST NOT BE SMALLER THAN X1 AND Y1')
            gv.error_Sound()
            sys.exit()
        if x2==x1 or y1==y2:
            print('-----------ERROR---------')
            print('X2 AND Y2 CANNOT EQUAL X1 AND Y1')
            gv.error_Sound()
            sys.exit()
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



    def do_Data_Analysis(self):
        #pretty broken right now
        self.catch_Errors_Analysis() #check if the user has any messed up input
        self.save_Settings()
        fileName=self.anlFileNameBox.get()
        folderPath=self.anlFolderPathBox.get()
        infoFile=open(folderPath+'\\'+fileName+'Info.txt')
        startVolt=float(infoFile.readline().split(',')[1])
        stopVolt=float(infoFile.readline().split(',')[1])
        infoFile.close()
        x1=int(self.x1BoxAnl.get())
        x2=int(self.x2BoxAnl.get())
        y1=int(self.y1BoxAnl.get())
        y2=int(self.y2BoxAnl.get())
        camName=self.cameraVarAnl.get() #which camera is being analyzed


        fitsFile = fits.open(folderPath+'\\'+fileName+camName+'.fits') #load the fits file into memory as a numpy array. This can
            #take a little while
        imagesArr = fitsFile[0].data #there is something called a header attached to a fit file. We only want the numpy arra
        fitsFile.close()
        imagesArr=np.flip(imagesArr,axis=1) #images are flipped when transferring from fits to numpy


        yMax=imagesArr.shape[1]
        y1N=yMax-y2
        if y1N<0:#to not loop back around if the user wants to use the whole image
            y1N=0
        y2N=yMax-y1
        imagesArr=imagesArr[:,y1N:y2N,x1:x2] #images is a 3 dimensional array where the first dimension is the images
            #the second is the rows (the y) and the third is the column (the x). Zero is the top left, not the bottom
            # right. It makes cropping a little tricky
        imagesSumArr=np.sum(np.sum(imagesArr, axis=2),axis=1) #sum along one axis and then the other. The result is an array where
            #each entry is the sum of all the pixels in that image.
        DAQData=np.loadtxt(folderPath+'\\'+fileName+'DAQData.csv',self.DAQDataArr,delimiter=',')
        MHzScaleArr=MakeMHzScale.make_MHz_Scale(DAQData)
        #MHzScaleArr=np.linspace(0,5000,num=DAQData.shape[0])
        galvoVoltArr=DAQData[:,0]
        liRefVoltArr=DAQData[:,1]
        P=np.polyfit(galvoVoltArr,MHzScaleArr, 1)
        temp=np.linspace(startVolt,stopVolt,num=imagesSumArr.shape[0])
        imageFreqMhzArr=P[1]+P[0]*temp


        #x1=MHzScaleArr
        #x2=imageFreqMhzArr

        #y1=liRefVoltArr
        #y1=-y1
        #y1=y1-np.min(y1)
        #y1=y1/np.max(y1)
        #y2=imagesSumArr
        #y2=y2-np.min(y2)
        #y2=y2/np.max(y2)
        #plt.plot(x1,y1)
        #plt.plot(x2,y2)
        #plt.show()

        S=float(self.saturationConstantBox.get()) #saturation constant, I/I_sat
        analyzer=Analyzer(imagesSumArr,imageFreqMhzArr,S)
        analyzer.fit_Image_Data()

        #print(analyzer.T)
#
#
        #self._Make_And_Save_Spectral_Fit_Plot(analyzer, DAQData)
        #plt.plot(imageFreqMhzArr,images)
        #plt.plot(imageFreqMhzArr,analyzer.fitFunc(imageFreqMhzArr))
        #plt.show()
    def _Make_And_Save_Spectral_Fit_Plot(self, analyzer, DAQData):
        fileName=self.anlFileNameBox.get()
        folderPath=self.anlFolderPathBox.get()
        galvoVoltArr=DAQData[:,0]
        liRefVoltArr=DAQData[:,1]
        plt.close('all')


        #now find the FWHM

        #xTemp=np.linspace(analyzer.imageFreqMHzArr[0],analyzer.imageFreqMHzArr[-1],num=10000)
        #yTemp=analyzer.fitFunc(xTemp)
        #yTemp=yTemp-yTemp.min()
        #yHalf=(yTemp.max()-yTemp.min())/2
        #xLeft=xTemp[np.argmax(yTemp>yHalf)]
        #xRight=xTemp[np.argmax(yTemp[np.argmax(yTemp):]<yHalf)+np.argmax(yTemp)]
        #FWHM=xRight-xLeft
        #print(xLeft,xRight)
#
        #plt.plot(xTemp,yTemp)
        #plt.show()





        x1=analyzer.imageFreqMHzArr
        y1=analyzer.imageAvgArr
        x1=np.linspace(x1[0],x1 [0],num)
        y2=analyzer.fitFunc(x) #the previously generated fit

        plt.figure( figsize=(10, 7))

        titleString="Temperature~ "+str(np.round(analyzer.T*1000))+' mk (dubious) '#+' FWHM: '+str(fwhm)+' mhz|'
        #titleString+='\n'+'Atom velocity: '+str(np.round(velocity, 1))+'m/s |'
        #titleString+=' signal size: '+str(sigSize)+' au | center frequency: '+str(np.round(F0, 1))+' mhz |'
        #titleString+='\n'+'region (x,y):'+str(region)+' | '+"Mhz per galvo volt= "+str(int(scale))+' |'+'\n'
        #titleString+=' created:'+str(datetime.date.today())

        plt.title(titleString)



        plt.ylabel("pixel value")
        plt.xlabel("MHz relative to F=2 transition")

        plt.plot(x,y1,label='Data')
        plt.plot(x,y2,label='Fit')
        plt.legend()
        plt.grid()
        #plt.savefig(folderPath+'\\'+fileName+self.cameraVarAnl.get())

        plt.show()



        #temp=round(1000*analyzer.T)
        #liRefForPlot=DAQData[:,1]
        ##create array of lithium PMT referance voltages to include with plot
        #liRefForPlot=liRefVoltage[start:stop+1]  #trim to only the value we care about
        #liRefForPlot=-liRefForPlot  #flip it cause it's negative
        #valueAdjust=(pixelValue.max()-pixelValue.min())/(liRefForPlot.max()-liRefForPlot.min())
        #liRefForPlot=.5*liRefForPlot*valueAdjust  #normalize to samae height as pixel value
        #liRefForPlot+=pixelValue.min()  #set the offset correctly
#
        #refFitForPlot=-refFit[start:stop+1]
        #refFitForPlot=.5*refFitForPlot*valueAdjust
        #refFitForPlot+=pixelValue.min()
#
        #analysisPlot=plt.figure(num=2, figsize=(10, 7))
        #plt.locator_params(nbins=20)
        #max=pixelValue.max()
        #ypos=(max-pixelValue.min())*.01+max
        #plt.axhline(y=max, color='black', linestyle='--')
        #plt.text(MHzScale[start], ypos, str(int(max)))
        #plt.axvline(x=F0, color='black', linestyle=':')  #centered on signal peak
        #plt.axvline(x=zeroPoint, color='black', linestyle=':', ymax=.5)  #centered on reference cell zero point
        #plt.plot(MHzScale[start:stop+1], pixelValue, label="pixel value")
        ##make high res fit
        #x=np.linspace(MHzScale[start], MHzScale[stop], num=10000)
        #y=dataAnalysis.spectralProfile(x, *PF)
        #plt.plot(x, y, label="fit")
        #plt.plot(MHzScale[start:stop+1], liRefForPlot, color='r', label="reference (AU)", alpha=.5)
        #plt.plot(MHzScale[start:stop+1], refFitForPlot, '--', color='g', label="reference fit", alpha=.5)
        #plt.xlabel("MHz relative to F=2 transition")
#
        #plt.ylabel("pixel value")
        #velocity=np.abs(gv.cLight*F0/(gv.Li_D2_Freq/1E6))
#
        #fitFunc=spi.interp1d(MHzScale[start:stop+1], pixelValue)  #to find fwhm
        #xNew=np.linspace(MHzScale[start], MHzScale[stop], 10000)
        #yNew=fitFunc(xNew)
        #yNew=yNew-np.average(yNew[:20])
        #yNew=-np.abs(yNew-(yNew.max()-yNew.min())/2)+(yNew.max()-yNew.min())*.01
        #peaks=sps.find_peaks(yNew, height=0)[0]
#
        #if peaks.shape[0]!=2:
        #    print("error! more than two peaks found during search for fwhm")
        #    plt.close()
        #    plt.plot(xNew, yNew)
        #    plt.show()
        #    sys.exit()
        #fwhm=np.round(np.abs(xNew[peaks[0]]-xNew[peaks[1]]), 1)
#
        #round_to_n=lambda x, n:round(x, -int(math.floor(math.log10(x)))+(n-1))
        #sigSize=round_to_n(pixelValue.max()-np.average(pixelValue[:10]), 3)
        #titleString=fitsFilePath+'\n'
        #titleString+="Temperature~ "+str(temp)+' mk (dubious) |'+' FWHM: '+str(fwhm)+' mhz|'
        #titleString+='\n'+'Atom velocity: '+str(np.round(velocity, 1))+'m/s |'
        #titleString+=' signal size: '+str(sigSize)+' au | center frequency: '+str(np.round(F0, 1))+' mhz |'
        #titleString+='\n'+'region (x,y):'+str(region)+' | '+"Mhz per galvo volt= "+str(int(scale))+' |'+'\n'
        #titleString+=' created:'+str(datetime.date.today())
        #if cheatVar.get()==1:
        #    titleString+="cheat method"
        #plt.title(titleString)
        #plt.tight_layout()
        #plotFileName=plotFilePath+"("+str(X1Anl_BOX.get())+','+str(X2Anl_BOX.get())+")"
        #plotFileName+="("+str(Y1Anl_BOX.get())+','+str(Y2Anl_BOX.get())+")"
        #plotFileName+=".png"
        #plt.grid()
        #plt.legend()
        #plt.savefig(plotFileName)
        #analysisPlot.show()

    def aquire_Data(self):
        self.catch_Errors_Aquisition()
        self.save_Settings()
        sweeper=Sweeper(self)
        sweeper.sweep()
        self.save_Fits_Files()
        self.save_Data()
        self.make_Info_File()
        #self.make_Lithium_Ref_Plot()

        #turn off the helium flow when done
        flowOut=DAQPin(gv.flowOutPin)
        flowOut.write(0)
        flowOut.close()
    def make_Lithium_Ref_Plot(self):
        fileName=self.dataFileNameBox.get()
        folderPath=self.dataFolderPathBox.get()
        x=self.DAQDataArr[:,0] #galvo voltage
        y=self.DAQDataArr[:,1] #lithium ref voltage
        plt.figure(figsize=(13,8))
        plt.plot(x,y)
        plt.title("Li Ref PMT Voltage vs Applied Galvo Voltage")
        plt.xlabel('Galvo Voltage, v')
        plt.ylabel('Lithium Ref PMT Voltage, v')
        plt.grid()
        plt.savefig(folderPath+'\\'+fileName+'LiRefVoltage')


    def make_Info_File(self):
        fileName=self.dataFileNameBox.get()
        folderPath=self.dataFolderPathBox.get()

        file=open(folderPath+'\\'+fileName+'Info.txt','w')
        file.write('Image aquisition start volt, '+str(self.startVoltBox.get())+'\n')
        file.write('Image aquisition stop volt, '+str(self.stopVoltBox.get())+'\n')
        file.write('Camera selection,'+str(self.cameraVarData.get())+'\n')
        file.write('Number of exposures(same for both cameras), '+str(self.expNumBox.get())+'\n')
        file.write('Exposure time Near, '+str(self.expTimeNearBox.get())+'\n')
        file.write('Exposure time Far, '+str(self.expTimeFarBox.get())+'\n')
        file.write('x1Near, '+str(self.x1NearBox.get())+'\n')
        file.write('x2Near, '+str(self.x2NearBox.get())+'\n')
        file.write('y1Near, '+str(self.y1NearBox.get())+'\n')
        file.write('y2Near, '+str(self.y2NearBox.get())+'\n')
        file.write('x1Far, '+str(self.x1FarBox.get())+'\n')
        file.write('x2Far, '+str(self.x2FarBox.get())+'\n')
        file.write('y1Far, '+str(self.y1FarBox.get())+'\n')
        file.write('y2Far, '+str(self.y2FarBox.get())+'\n')
        file.write('Data and time of file creation (year,month,day,time), '+str(datetime.now())+'\n')

        file.close()


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

    def make_Flow(self):
        b=4.29
        m=1.4236
        flowOut=DAQPin(gv.flowOutPin)
        flowDesired=float(self.flowRateBox.get())
        if flowDesired>0.0:
            flowOut.write((5.0/2000)*(flowDesired-b)/m)
        else:
            flowOut.write(0.0)
        flowOut.close(zero=False)
    def Li_Reference_Check(self):
        plt.close('all')
        voltArr=np.linspace(gv.minScanVal,gv.maxScanVal,num=250)
        temp=[]
        galvoOut=DAQPin(gv.galvoOutPin)
        liRefIn=DAQPin(gv.lithiumRefInPin)
        for volt in voltArr:
            galvoOut.write(volt)
            temp.append(liRefIn.read(numSamples=1000))
        galvoOut.close()
        liRefIn.close()
        plt.plot(voltArr,temp)
        plt.grid()
        plt.title('Lithium reference chamber signal vs galvo voltage')
        plt.xlabel('Galvo voltage, v')
        plt.ylabel('PMT voltage, v')
        plt.show()
    def load_Settings(self):
        try:
            file=open("GUI_Settings.txt", "r")
        except:
            print("NO SETTINGS FILE FOUND")
            return
        i=0
        for item in file:
            item=item.strip()
            if i>=len(self.settingsList):
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