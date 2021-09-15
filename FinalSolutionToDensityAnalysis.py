from UncertaintyAnalysis_Functions import find_Chi_Squared_Uncertainty
from DataAnalysis import fit_Spectral_Data
import globalVariables as gv
import numpy as np
import os
import matplotlib.pyplot as plt
import scipy.special as sps
from scipy.integrate import quad
import scipy.optimize as spo
from phaseSpaceAnalysis_Functions import get_Images_Array_And_MHz_Scale_Arr
from phaseSpaceAnalysis_Functions import subtract_Background
from phaseSpaceAnalysis_Functions import voigtImageFit

#Computes the peak density of a Fluorescence Run. Requires user to enter many parameters!
#ALL values are in SI unless specifically stated otherwise

#Directory
path = "C:\Data\Runs\\8_29_21"
os.chdir(path)
fileName = 'run30Far'

imagesArr,imageFreqMhzArr = get_Images_Array_And_MHz_Scale_Arr(path,fileName)
imagesArrOriginal = imagesArr
imagesArr = subtract_Background(imagesArr,3)


# testArr = np.sum(imagesArr[:,0:-1,0:-1],axis=0)
#
# for x in range (0,np.shape(testArr)[1]):
#     max = np.max(testArr[110:116,x])
#     for y in range (0,np.shape(testArr)[0]):
#         value = testArr[y,x]
#         testArr[y,x] = value/max
#
# plt.imshow(testArr[70:150,0:-1],vmin = 0.0, vmax=1.0)
# plt.show()


import time
t=time.time()

showPlots = True
lensHeating = True
#---Parameters Relating to Optics/Camera Setup---

#Laser Parameters
total_laser_power = 243e-9 #Units of W
wx = 5.900e-4 #Beam waist of laser in transverse direction in m
wz = 2.842e-2 #Beam waist of laser in longitudinal direction in m
laser_profile_center = 98.5 #Location of peak laser intensity in terms of pixels.
laserJitter=1.1 #MHz
window_efficiency = 0.93
polarization_correction_factor = 2

#Camera Parameters
exposure_time = 1 #seconds
binning = 4
pixel_size_on_camera = 24e-6 #m
magnification = 3.0 #Measured by placing ruler in front of camera
quantum_efficiency = 0.6
lens_diameter = 0.048 #m
lens_to_laser_length = 0.152 #m
counts_to_e = 0.5

AtomVelocity = 20700 #cm/s
pixelSize = magnification*binning*pixel_size_on_camera


#Size of nxn pixel box
boxSize = 2

xstart = 50
xdelta = 140

ystart = 90
ydelta = 40

if showPlots:
    image=np.mean(imagesArr[0:-1, 0:-1, 0:-1], axis=0)
    plt.imshow(image)
    plt.title('Average of Pixel Values over All Images')
    plt.xlabel('Longitudinal Direction [Pixel Number]')
    plt.ylabel('Vertical Direction [Pixel Number]')
    plt.axhline(y=ystart,c='r',linestyle="--")
    plt.axhline(y=ystart+ydelta,c='r',linestyle="--")
    plt.axvline(x=xstart,c='r',linestyle="--")
    plt.axvline(x=xstart+xdelta,c='r',linestyle="--")
    plt.show()

xpixels = []
ypixels = []
images = []
total_counts = []

for image in range (0,len(imagesArr)):

    for z in range (xstart,xstart+xdelta,boxSize):

        for y in range (ystart,ystart+ydelta,boxSize):
            temp = np.sum(imagesArr[image,y:y+boxSize,z:z+boxSize])
            xpixels.append(z)
            ypixels.append(y)
            images.append(image)
            total_counts.append(temp)

maxElement = total_counts.index(max(total_counts))
xMaxTemp = xpixels[maxElement]
yMax = ypixels[maxElement]
imageMax = images[maxElement]
countsMax = total_counts[maxElement]

#Test what max pixel signal is nearby xMax. Should be below 60,000 to not saturate camera.
TestSize = 20
max_pixel_value = np.max((imagesArr[imageMax,yMax-TestSize:yMax+TestSize,xMaxTemp-TestSize:xMaxTemp+TestSize]))


def bigAssFunction(xMax,countsMax,truth):

    #----Frequency Fit Over Entire Height----

    freqSignal=np.sum(np.sum(imagesArr[:, 0:-1, xMax:xMax+boxSize], axis=2), axis=1)
    x=np.arange(freqSignal.shape[0])

    #----Computing Peak Density----
    #Computing S(w)
    area_under_curve=np.trapz(freqSignal, x=2*np.pi*imageFreqMhzArr*10**6)
    maxSignal=np.max(freqSignal)
    S_w=maxSignal/area_under_curve  #s/rad

    #Compute frequency integrated cross section
    cross_section=(3/4)*gv.lambda_D2line**2*gv.A21*S_w

    #Compute fraction of light collected and total number of photons per second
    fractional_solid_angle=(lens_diameter/2)**2/(4*lens_to_laser_length**2)
    photons_per_second=countsMax*(1/exposure_time)*(1/quantum_efficiency)*(1/counts_to_e)\
                       *(1/window_efficiency)*(1/fractional_solid_angle)

    #Integration of laser intensity over volume
    IntegrationConst=2*total_laser_power*window_efficiency/(np.pi*wx*wz)
    yIntegration=boxSize*pixelSize
    xIntegration=np.sqrt(np.pi/2)*wx

    def integrand(z):
        return np.exp(-2*z**2/wz**2)

    z1Bound=xMax*pixelSize-laser_profile_center*pixelSize
    z2Bound=(xMax+boxSize)*pixelSize-laser_profile_center*pixelSize
    zIntegration=quad(integrand, z1Bound, z2Bound)
    I=IntegrationConst*yIntegration*xIntegration*zIntegration[0]

    peak_density=gv.hbar*gv.w21*photons_per_second*polarization_correction_factor/(I*cross_section)/1e6  #cm^-3

    #----Spatial Profile----

    def voigtSpaceFit(x, x0, a, sigma, gamma, b):
        v0=sps.voigt_profile(0, sigma, gamma)
        v=sps.voigt_profile(x-x0, sigma, gamma)
        return a*v/v0+b

    #This is in mm just because
    spatial_signal=np.sum(np.sum(imagesArr[:, 0:-1, xMax:xMax+boxSize], axis=0), axis=1)
    spatial_y=np.arange(0, spatial_signal.shape[0])*pixelSize*1000

    eps=1e-10  #small number
    guess=[spatial_y[np.argmax(spatial_signal)], spatial_signal.max()-spatial_signal.min(), 1.0, 1.0,
           spatial_signal.min()]

    bounds=[(-np.inf, eps, eps, eps, -np.inf), (np.inf, np.inf, np.inf, np.inf, np.inf)]

    spatial_params=spo.curve_fit(voigtSpaceFit, spatial_y, spatial_signal, p0=guess, bounds=bounds)[0]

    x0, a, sigma, gamma, b=spatial_params  #gamma is HWHM NOT FWHM

    spatial_fwhm=.5346*(2*spatial_params[3])+np.sqrt(.2166*(2*spatial_params[3])**2+(spatial_params[2]*2.335)**2)

    #----Peak Intensity and Flux in 1cm Circle----
    peak_intensity=peak_density*AtomVelocity  #units of cm^-2s^-1

    spatial_gamma=2*gamma/10  #this is the FWHM in cm NOT HWHM in mm

    def integrand(z):
        return z*(spatial_gamma/2)**2/(z**2+(spatial_gamma/2)**2)

    radius=0.5  #cm
    flux=peak_intensity*2*np.pi*quad(integrand, 0, radius)[0]

    if truth == True:

        print()
        print('|---Image and Box Location for Max Signal---|')
        print('|Max Pixel Value:', max_pixel_value)
        print('|Image on Fits File with Max Signal:', imageMax+1)
        print('|x pixel coordinates:', xMax, xMax+boxSize)
        print('|y pixel coordinates:', yMax, yMax+boxSize)
        print('|Max Total Counts:', countsMax)
        print()

        print('Spatial gamma, FWHM (mm):',2*spatial_params[3])

        if showPlots:
            plt.scatter(spatial_y, spatial_signal, c='r', label='data')
            plt.plot(spatial_y, voigtSpaceFit(spatial_y, *spatial_params), label='fit')
            plt.suptitle('Spatial Profile')
            #plt.title('Voigt FWHM ='+str(np.round(spatial_fwhm, 2))+'mm')
            plt.ylabel('Intensity (a.u.)')
            plt.xlabel('Transverse Position (mm)')
            plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
            #plt.grid()
            plt.legend()
            plt.show()

        #---Uncertainty in Spatial Gamma---
        optimal_gamma = spatial_params[3]
        test_size = 0.4
        chi_squared = []
        test_gamma = []
        for gamma_value in np.arange(optimal_gamma-test_size,optimal_gamma+test_size,test_size/50):

            test_gamma.append(gamma_value)
            def spatialVoigtGammaFixed(x, x0, a, sigma, b):
                gamma=gamma_value
                return voigtSpaceFit(x,x0,a,sigma,gamma,b)

            guess_fixed=[spatial_y[np.argmax(spatial_signal)], spatial_signal.max()-spatial_signal.min(), 1.0,
                   spatial_signal.min()]

            bounds_fixed=[(-np.inf, eps, eps, -np.inf), (np.inf, np.inf, np.inf, np.inf)]

            spatial_params_fixed=spo.curve_fit(spatialVoigtGammaFixed, spatial_y, spatial_signal,
                                               p0=guess_fixed,bounds = bounds_fixed)[0]

            R2 = []
            for i in np.arange(0,len(spatial_y)):
                value_fit=spatialVoigtGammaFixed(spatial_y[i],*spatial_params_fixed)
                value_data = spatial_signal[i]
                difference = (value_fit-value_data)**2
                R2.append(difference)

            chi_squared.append(np.sum(R2))
            # print('chi squared:',chi_squared)
            # print(gamma_value)
            #
            # plt.scatter(spatial_y, spatial_signal, c='r', label='data')
            # plt.plot(spatial_y, spatialVoigtGammaFixed(spatial_y, *spatial_params_fixed), label='fit')
            # plt.suptitle('Spatial Profile')
            # #plt.title('Voigt FWHM ='+str(np.round(spatial_fwhm, 2))+'mm')
            # plt.ylabel('Intensity (a.u.)')
            # plt.xlabel('Transverse Position (mm)')
            # plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
            #
            # plt.legend()
            # plt.show()

        find_Chi_Squared_Uncertainty(test_gamma,chi_squared)
        #---End of Spatial Uncertainty---

        #---Uncertainty in S(w)---
        backgroundSignal = []
        for I in range(0,3):
            test= np.sum(imagesArrOriginal[I,10:11,xMax:xMax+boxSize])
            backgroundSignal.append(test)
            print('image',test)

        STDEV = np.std(backgroundSignal)
        stdev_mean = STDEV / np.sqrt(len(backgroundSignal))*256
        print(stdev_mean)

        #----Frequency Fit Over Entire Height----

        freqSignal=np.sum(np.sum(imagesArr[:, 0:-1, xMax:xMax+boxSize], axis=2), axis=1)+stdev_mean
        x=np.arange(freqSignal.shape[0])

            #----Computing Peak Density----
        #Computing S(w)
        area_under_curve=np.trapz(freqSignal, x=2*np.pi*imageFreqMhzArr*10**6*1.1)
        maxSignal=np.max(freqSignal)
        S_w_test=maxSignal/area_under_curve  #s/rad

        #Compute frequency integrated cross section
        cross_section_test=(3/4)*gv.lambda_D2line**2*gv.A21*S_w_test
        print('S(w) upper bound:',S_w_test)
        print('Cross section upper bound:',cross_section_test)

        #---End of S(w) Uncertainty---
        guess=[np.argmax(freqSignal), np.max(freqSignal)-np.min(freqSignal),1, 1, np.min(freqSignal)]
        params, pcov=spo.curve_fit(voigtImageFit, x, freqSignal, p0=guess)

        spectralFit=fit_Spectral_Data(imageFreqMhzArr, freqSignal, lensHeating=lensHeating, vTMaxLens=.07*200*.9
                                      , peakMode='single', laserJitter=laserJitter)
        spectralFit.print_Results()

        spectralFit=fit_Spectral_Data(imageFreqMhzArr, freqSignal, lensHeating=False, vTMaxLens=.07*200*.9
                                      , peakMode='single', laserJitter=laserJitter)
        spectralFit.print_Results()

        if showPlots:
            plt.scatter(imageFreqMhzArr, freqSignal, c='r',label='Data')
            plt.plot(imageFreqMhzArr,voigtImageFit(x, *params),label='Fit')
            plt.suptitle('Summed Signal vs. Image')
            plt.xlabel('Relative Frequency (MHz)')
            plt.ylabel('Intensity (a.u.)')
            plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
            plt.legend()
            plt.show()

        print()
        print('|-----Summary of Results-----|')
        print('|S(w0) [s/rad]', maxSignal/area_under_curve)
        print('|Photons per second:', "{:e}".format(photons_per_second))
        print('|Laser Integration [Wm]:', "{:e}".format(I))
        print('|Cross Section [m^2]', "{:e}".format(cross_section))
        print('|Spatial Voigt FWHM [mm]:', spatial_fwhm)
        print('|Peak Density [cm^-3]:', "{:e}".format(peak_density))
        print('|Peak Intensity [cm^-2s^-1]:',"{:e}".format(peak_intensity))
        print('|Flux in 1cm circle [atoms/s]:', "{:e}".format(flux))

        if max_pixel_value>=50000:
            print('WARNING!! CAMERA IS SATURATED!!')
            print('MAX PIXEL VALUE:',max_pixel_value)

        print('t:', time.time()-t)

    else :
        return flux

fluxArr = []
xPositions = []
counts = []
for x in range (xstart,xstart+xdelta,boxSize):

    total_counts = np.sum(imagesArr[imageMax,yMax:yMax+boxSize,x:x+boxSize])
    counts.append(total_counts)
    xPositions.append(x)
    fluxArr.append(bigAssFunction(x,total_counts,False))

plt.scatter(xPositions,fluxArr)
plt.xlabel('Longitudinal Position [pixel #]')
plt.ylabel('Flux [atoms/s]')
plt.title('Flux vs. Longitudinal Position')
plt.show()
maxElement = fluxArr.index(max(fluxArr))
xMax = xPositions[maxElement]
countsMax = counts[maxElement]

bigAssFunction(xMax,countsMax,True)