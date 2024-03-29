import time
import globalVariables as gv
from DAQClass import DAQPin
import numpy as np
import matplotlib.pyplot as plt

#Compute the spectroscopy laser linewidth using the 1 GHz FSR Etalon.
#Method is to measure the voltage between the peaks to get a MHz per voltage. Next measure the slope at the approx. linear
#portion of the curve. Then park the laser at this point and measure the rms signal over some time scale.
#The rms of the voltage signal at the linear portion can then be converted to a MHz rms.

#Adjustable values. Datapoints is the number of points to construct the etalon's signal. Samples is the number of times
#the RMS signal is collected with pauses inbetween. Offset is the number of datapoints to the left/right of the
#half max of of the peak. Make sure that the offset isn't so large that the linear approx. breaks down.
datapoints = 100
samples = 10
offset = 2
wait_time = 1 #wait time in (s) between samples

#FSR of the Elaton
FSR = 1000

#Construct the Etalon's signal over a designated voltage range.
galvoOut=DAQPin(gv.galvoOutPin)
voltArr=np.linspace(-0.75,1.25 ,num = datapoints)
signalArr = []
laserPin=DAQPin(gv.laserWidthInPin)
for volt in voltArr:
    galvoOut.write(volt)
    time.sleep(0.01)
    signalArr.append(laserPin.read(numSamples=1000,average=True))
    print(laserPin.read(numSamples=1000,average=True))
laserPin.close()
galvoOut.close()


#Find the max of each of the two peaks in terms of x and y voltages as well as the min value between them.
midway = int(datapoints / 2)
FirstPeak = np.argmax(signalArr[0:midway])
SecondPeak = midway + np.argmax(signalArr[midway:datapoints])
MinBetweenPeaks = FirstPeak + np.argmin(signalArr[FirstPeak:SecondPeak])
DeltaPeaks = SecondPeak - FirstPeak
FirstPeakVoltage = signalArr[FirstPeak]
SecondPeakVoltage = signalArr[SecondPeak]
MidwayVoltage = signalArr[MinBetweenPeaks]
LinearSlopeCenterVoltage = (FirstPeakVoltage + MidwayVoltage) / 2
print('LinearSlopeCenterVoltage',LinearSlopeCenterVoltage)


#Contruct the Mhz per voltage scale.
MHzScale = FSR / (voltArr[SecondPeak]-voltArr[FirstPeak])

def find_nearest(array, value):
    array=np.asarray(array)
    idx=(np.abs(array-value)).argmin()
    return array[idx]

def binarySearchCount(arr, n, key):
    left=0
    right=n

    mid=0
    while (left<right):
        mid=(right+left)//2

        if (arr[mid]==key):

            while (mid+1<n and arr[mid+1]==key):
                mid+=1
            break

        elif (arr[mid]>key):
            right=mid

        else:
            left=mid+1

    while (mid>-1 and arr[mid]>key):
        mid-=1

    return mid+1


#Determine the center location of the linear portion of the peak and the(x1,y1) and (x2,y2) values to find the slope.
LinearSlopeCenterVoltage = (FirstPeakVoltage + MidwayVoltage) / 2
LinearSlopeLocation = find_nearest(signalArr[FirstPeak:MinBetweenPeaks],LinearSlopeCenterVoltage)
LinearSlopeIndex = signalArr.index(LinearSlopeLocation)
deltaV = voltArr[LinearSlopeIndex+offset] - voltArr[LinearSlopeIndex-offset]
Slope = deltaV * MHzScale / (signalArr[LinearSlopeIndex-offset] - signalArr[LinearSlopeIndex+offset])


#Measure the rms signal at the center point of the linear portion x number of times (x=samples) and average the results
#as well as compute the SD. This should give an estimate as to how the laser linewidth changes over long time scales.
galvoOut=DAQPin(gv.galvoOutPin)
laserPin=DAQPin(gv.laserWidthInPin)
galvoOut.write(voltArr[LinearSlopeIndex])

RMS =[]
for i in range(samples):
    RMS.append(laserPin.read(numSamples=250000,average=False,error = True))
    time.sleep(wait_time)


Noise = laserPin.read(numSamples=250000,average=4,error=3)
folderPath='C:\Data\Runs\8_23_21\\'
np.savetxt(folderPath+'LaserNoiseFreq ',(Noise*MHzScale))
Order = np.sort(Noise)
NBins = 75
Bins = np.linspace(Order[0],Order[-1],num= NBins)



NoiseValues = []
i = 0

while i < NBins:
    if i==0:
        value =  binarySearchCount(Order,len(Order),Bins[i])
    else :
        value = binarySearchCount(Order,len(Order),Bins[i])-binarySearchCount(Order,len(Order),Bins[i-1])
    NoiseValues.append(value)
    i = i + 1


RMS_average = np.average(RMS)
RMS_SD = np.std(RMS)


#Measure the signal with the laser set to the mininum between the peaks. Laser Fluctuations should have minimal
#effect on the etalon's signal so the fluctuations are primarily noise from the Etalon. This was tested by blocking
#the laser input to the etalon and rms of the signals is within a few percent.
galvoOut.write(voltArr[MinBetweenPeaks])
DAQNoise = laserPin.read(numSamples=250000,average=False,error=True)
laserPin.close()
galvoOut.close()



#Compute the RMS linewidth and then adjusted linewidth assuming the laser linewidth and DAQ board noise adds in quadrature.
Linewidth = Slope * RMS_average
AdjustedLinewidth = np.sqrt(RMS_average**2-DAQNoise**2) * Slope
AdjustedLinewidth_SD = np.sqrt((RMS_average+RMS_SD)**2-DAQNoise**2) * Slope - AdjustedLinewidth

print('MHz per Volt', MHzScale)
print('delta x',deltaV)
print('delta y', signalArr[LinearSlopeIndex-offset] - signalArr[LinearSlopeIndex+offset])
print('RMS', RMS_average)
print('RMS SD', RMS_SD)
print('DAQ Noise', DAQNoise)
print('Linewidth (MHz):',np.round(Linewidth,2))
print('Linewidth adjusted for DAQ Noise (MHz):',np.round(AdjustedLinewidth,2))
print('Adjusted Linewidth SD (MHz)', np.round(AdjustedLinewidth_SD,3))

plt.plot(voltArr,signalArr)
plt.axvline(voltArr[FirstPeak])
plt.axvline(voltArr[SecondPeak])
plt.axvline(voltArr[LinearSlopeIndex-offset])
plt.axvline(voltArr[LinearSlopeIndex+offset])
plt.show()

plt.scatter(Bins,NoiseValues)
plt.title('Laser Noise Profile')
plt.show()

