import sys
import numpy as np
import tkinter as tk
import matplotlib.pyplot as plt
import os


class GUI:
    def __init__(self):
        self.settingsList=[]
        self.x1=None
        self.x2=None
        self.y1=None
        self.y2=None
        self.voltArr=None  #array to to hold voltage value sto scan over
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

        self.numImages=tk.Entry(self.window)
        self.numImages.config(width=5)
        self.numImages.grid(column=1, row=1, sticky='W')
        self.settingsList.append(self.numImages)

        lbl31=tk.Label(self.window, text='Exp time (ms)')
        lbl31.grid(column=2, row=1)

        self.expTime=tk.Entry(self.window)
        self.expTime.config(width=5)
        self.expTime.grid(column=3, row=1, sticky='W')
        self.settingsList.append(self.expTime)

        lbl32=tk.Label(self.window, text='bin size')
        lbl32.grid(column=4, row=1)

        self.binSize=tk.Entry(self.window)
        self.binSize.config(width=5)
        self.binSize.grid(column=5, row=1, sticky='W')
        self.settingsList.append(self.binSize)

        lbl4=tk.Label(self.window, text='Camera')
        lbl4.grid(column=0, row=4)

        self.cameraVarData=tk.StringVar(self.window)
        self.cameraVarData.set("FAR")
        cameraChoice=["NEAR", "FAR"]
        CAMERA_MENU=tk.OptionMenu(self.window, self.cameraVarData, *cameraChoice)
        CAMERA_MENU.grid(column=1, row=4, columnspan=2)
        self.settingsList.append(self.cameraVarData)

        lbl5=tk.Label(self.window, text='image x1')
        lbl5.grid(column=0, row=5)
        self.x1=tk.Entry(self.window)
        self.x1.config(width=5)
        self.x1.grid(column=1, row=5, sticky='W')
        self.settingsList.append(self.x1)

        lbl5=tk.Label(self.window, text='image x2')
        lbl5.grid(column=2, row=5)
        self.x2=tk.Entry(self.window)
        self.x2.config(width=5)
        self.x2.grid(column=3, row=5, sticky='W')
        self.settingsList.append(self.x2)

        lbl5=tk.Label(self.window, text='image y1')
        lbl5.grid(column=0, row=6)
        self.y1=tk.Entry(self.window)
        self.y1.config(width=5)
        self.y1.grid(column=1, row=6, sticky='W')
        self.settingsList.append(self.y1)

        lbl5=tk.Label(self.window, text='image y2')
        lbl5.grid(column=2, row=6)
        self.y2=tk.Entry(self.window)
        self.y2.config(width=5)
        self.y2.grid(column=3, row=6, sticky='W')
        self.settingsList.append(self.y2)

        self.saveDataVar=tk.BooleanVar()
        saveDataCheckButton=tk.Checkbutton(self.window, text='save data', variable=self.saveDataVar)
        saveDataCheckButton.grid(column=1, row=7)
        self.settingsList.append(self.saveDataVar)

        self.chopperVar=tk.BooleanVar()
        chopperVarButton=tk.Checkbutton(self.window, text='chopper', variable=self.chopperVar)
        chopperVarButton.grid(column=1, row=8)
        self.settingsList.append(self.chopperVar)

        self.showPlotVar=tk.BooleanVar()
        showDataAnalysiButton=tk.Checkbutton(self.window, text='Show plot',
                                             variable=self.showPlotVar)
        showDataAnalysiButton.grid(column=1, row=9)
        self.settingsList.append(self.showPlotVar)

        lbl3=tk.Label(self.window, text='Run name')
        lbl3.grid(column=0, row=15)

        self.fileName=tk.Entry(self.window)
        self.fileName.config(width=5)
        self.fileName.grid(column=1, row=15, sticky='W')
        self.settingsList.append(self.fileName)

        lbl3=tk.Label(self.window, text='Folder path')
        lbl3.grid(column=0, row=16)

        self.folderPath=tk.Entry(self.window)
        self.folderPath.config(width=5)
        self.folderPath.grid(column=1, row=16, sticky='W')
        self.settingsList.append(self.folderPath)

        runButton=tk.Button(self.window, text='RUN', font=("Arial", 20), background="green", command=self.run)
        runButton.config(height=2, width=10)
        runButton.grid(column=0, row=17, columnspan=4, rowspan=4)
        self.load_Settings()

        self.window.protocol("WM_DELETE_WINDOW", self.close_GUI)
        self.window.mainloop()

    def close_GUI(self):
        self.save_Settings()
        self.window.destroy()
        sys.exit()

    def run(self):
        self.save_Settings()
        self.voltArr=np.linspace(int(self.voltStartBox.get()), int(self.voltStopBox.get()),
                                 num=int(self.numImages.get()))
        if os.path.isdir(self.folderPath.get())==False:
            print('-----ERROR-----------')
            print('YOU HAVE ENTERED AN INVALID FOLDERPATH')
        if self.chopperVar.get()==True:
            self.run_With_Chopper()

    def take_Dark_Image_Average(self, num=3):
        None

    def run_With_Chopper(self):
        y1=(np.pi)*np.exp(-self.voltArr**2)  #chopper on
        y2=np.exp(-self.voltArr**2)  #chopper off

        y1Int=np.trapz(y1)
        y2Int=np.trapz(y2)
        ratio=y1Int/y2Int

        plt.suptitle('Signal with chopper on (closed) and off (open)')
        plt.title('ratio of integral of on to off = '+str(np.round(ratio, 3)))
        plt.plot(self.voltArr, y1, label='Chopper on')
        plt.plot(self.voltArr, y2, label='Chopper off')
        plt.xlabel('Volts')
        plt.ylabel('Pixel counts')
        plt.legend()
        plt.grid()
        if self.saveDataVar.get()==True:
            #save file
            None
        if self.showPlotVar.get()==True:
            plt.show()

    def run_Without_Chopper(self):
        y=np.exp(-self.voltArr**2)
        if self.saveDataVar.get()==True:
            #save file
            None

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