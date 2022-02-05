from UncertaintyAnalysis_Functions import find_Chi_Squared_Uncertainty
from UncertaintyAnalysis_Functions import FluorDensityUncertainty
from UncertaintyAnalysis_Functions import f_uncertainty
from UncertaintyAnalysis_Functions import I_uncertainty
from UncertaintyAnalysis_Functions import z_uncertainty
from UncertaintyAnalysis_Functions import delta_flux
from UncertaintyAnalysis_Functions import voigtFWHM_uncertainty
from UncertaintyAnalysis_Functions import brightness_uncertainty
from UncertaintyAnalysis_Functions import qGaussFWHM_uncertainty
from DataAnalysis import fit_Spectral_Data
from NewFile import fit_Spectral_Data_ConstrainedSigma
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
from qGaussian import qGaussianHelper

#Computes the peak density of a Fluorescence Run. Requires user to enter many parameters!
#ALL values are in SI unless specifically stated otherwise

#Directory
path = "C:\Data\Runs\\1_27_22"
os.chdir(path)
fileName = 'run56Far'

imagesArr,imageFreqMhzArr = get_Images_Array_And_MHz_Scale_Arr(path,fileName)
print('first and last freq', imageFreqMhzArr[0], imageFreqMhzArr[-1])
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

#Camera Parameters
exposure_time = 1 #seconds
binning = 4
pixel_size_on_camera = 24e-6 #m
magnification = 3.0 #Measured by placing ruler in front of camera
M_uncertainty = 0.02
quantum_efficiency = 0.6
lens_diameter = 0.048 #m
diameter_uncertainty = 0.0005 #m
lens_to_laser_length = 0.152 #m
laser_length_uncertainty = 0.003 #m
counts_to_e = 0.5

#Laser Parameters
total_laser_power = 200e-9 #Units of W
power_uncertainty = total_laser_power * 0.07
wx = 5.900e-4 #Beam waist of laser in transverse direction in m
wz = 2.842e-2 #Beam waist of laser in longitudinal direction in m
wz_uncertainty = 0.0006 #uncertainty in m
laser_profile_center = 107.4 * binning / 4 #Location of peak laser intensity in terms of pixels.
laser_center_uncertainty = 2.4 * binning / 4 #uncertainty in pixels
laserJitter=1.1 #MHz
window_efficiency = 0.93
eff_uncertainty = window_efficiency * 0.025
polarization_correction_factor = 2

AtomVelocity = 21000 #cm/s
AtomVelocity_uncertainty = 200 #cm/s
pixelSize = magnification*binning*pixel_size_on_camera #corrected pixel size. aka "Superpixel" size.


#Size of nxn pixel box
boxSize = 4

xstart = 40 #originally 90
xdelta = 140

ystart = 70
ydelta = 80

assert imagesArr[0].shape == (int(256 * 4/ binning),int(256 * 4/ binning))

if showPlots:
    image=np.mean(imagesArr[0:-1, 0:-1, 0:-1], axis=0)
    plt.imshow(image)
    plt.title('Average of Pixel Values over All Images')
    plt.xlabel('Longitudinal Direction (Pixel Number)',fontsize=11)
    plt.ylabel('Transverse Direction (Pixel Number)',fontsize=11)
    plt.axhline(y=ystart,c='r',linestyle="--")
    plt.axhline(y=ystart+ydelta,c='r',linestyle="--")
    plt.axvline(x=xstart,c='r',linestyle="--")
    plt.axvline(x=xstart+xdelta,c='r',linestyle="--")
    plt.colorbar()
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

    y0 = boxSize*pixelSize
    z1Bound=(xMax-laser_profile_center)*pixelSize

    z2Bound=((xMax+boxSize)-laser_profile_center)*pixelSize

    I = 0.5*total_laser_power*y0*(sps.erf(np.sqrt(2)*z2Bound / wz)-sps.erf(np.sqrt(2)*z1Bound / wz))

    const = gv.hbar*gv.w21 * polarization_correction_factor / (exposure_time * quantum_efficiency * counts_to_e)
    f = 16 * lens_to_laser_length**2 / lens_diameter**2
    C = countsMax
    eff = window_efficiency

    peak_density = const * C * f / (eff**2 * I * cross_section) / 1e6

    #----Spatial Profile----

    def voigtFit(x, x0, a, sigma, gamma, b):
        v0=sps.voigt_profile(0, sigma, gamma)
        v=sps.voigt_profile(x-x0, sigma, gamma)
        return a*v/v0+b

    #---For Angled Beam DATA ONLY
    # guess_angled=[0, np.max(freqSignal)-np.min(freqSignal), 1, 1, np.min(freqSignal)]
    # angledParams, pcov=spo.curve_fit(voigtFit, imageFreqMhzArr, freqSignal, p0=guess_angled)
    #
    # angledGamma=angledParams[3]
    # angledSigma=angledParams[2]
    # print(angledGamma, angledSigma)
    # angledSpectralFWHM=0.5346*2*angledGamma+np.sqrt(0.2166*4*angledGamma**2+2.35**2*angledSigma**2)
    # print(angledSpectralFWHM)
    # V_longitudinal_spread=gv.cLight*angledSpectralFWHM*10**6*2*np.pi/gv.w21
    # print(V_longitudinal_spread)
    # plt.scatter(imageFreqMhzArr, freqSignal)
    # plt.plot(imageFreqMhzArr, voigtFit(imageFreqMhzArr, *angledParams))
    # plt.show()

    #---END OF ANGLED BEAM

    #This is in units of mm
    spatial_signal=np.sum(np.sum(imagesArr[:, 0:-1, xMax:xMax+boxSize], axis=0), axis=1)
    spatial_signal = spatial_signal / np.max(spatial_signal)#normalize height to 1
    #spatial_signal= (spatial_signal/ 1508515.6666666672)  * (153/180.7) # normalize to other runs

    spatial_y=np.arange(0, spatial_signal.shape[0])*pixelSize*1000

    #Using Voigt Profile
    eps=1e-10  #small number
    guess=[spatial_y[np.argmax(spatial_signal)], spatial_signal.max()-spatial_signal.min(), 1.0, 1.0,
           spatial_signal.min()]
    bounds=[(-np.inf, eps, eps, eps, -np.inf), (np.inf, np.inf, np.inf, np.inf, np.inf)]
    spatial_params=spo.curve_fit(voigtFit, spatial_y, spatial_signal, p0=guess, bounds=bounds)[0]
    x0, a, sigma, gamma, b=spatial_params  #gamma is HWHM NOT FWHM
    #order of parameters are switched for voigt function definied internally and that in data analysis. Should fix!
    spatial_FWHM=.5346*(2*spatial_params[3])+np.sqrt(.2166*(2*spatial_params[3])**2+(spatial_params[2]*2.335)**2)


    #---Using qGaussian Profile---
    qG=qGaussianHelper()

    qGguessParams=[spatial_y[np.argmax(spatial_signal)], spatial_signal.max()-spatial_signal.min(),
                   spatial_signal.min(), 0.25, 2, 1.0]

    qGbounds=[(-np.inf, 1e-9, -np.inf, 1e-6, 1e-9, 1e-9), (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf)]
    qGparams, pcov=spo.curve_fit(qG, spatial_y, spatial_signal, bounds= qGbounds, p0=qGguessParams)

    # plt.scatter(spatial_y, spatial_signal, c = 'r')
    # plt.plot(spatial_y, qG(spatial_y, *qGparams))
    # plt.title('Spatial profile fit to q Guassian ')
    # plt.show()

    spatialFWHM = qG.get_FWHM(*qGparams)
    # print(spatialFWHM)

    qG_sigma = qGparams[3]

    radius = 5 #units of mm
    # spatial_x_Values = np.linspace(-3*qG_sigma, 3*qG_sigma, 10_000)

    left_x = np.linspace(-radius, 0, 10_000)
    right_x = np.linspace(0, radius, 10_0000)

    left_FitValues = qG(left_x,0,1.0,0, qG_sigma, qGparams[4], qGparams[5])
    right_FitValues = qG(right_x,0,1.0,0, qG_sigma, qGparams[4], qGparams[5])

    left_integrand = np.abs(np.multiply(left_x, left_FitValues))
    right_integrand = np.abs(np.multiply(right_x, right_FitValues))
    left_flux_integral = np.trapz(left_integrand, left_x)
    right_flux_integral = np.trapz(right_integrand, right_x)

    flux_integral = np.pi*(left_flux_integral + right_flux_integral)/100 #factor of 100 to convert mm^2 to cm^2

    # spatial_params = spo.curve_fit(qGaussian, spatial_y, spatial_signal, p0=guess)[0]
    #
    # plt.scatter(spatial_y, spatial_signal, c = 'r')
    # plt.plot(spatial_y, qGaussian(spatial_y, *spatial_params))
    # plt.show()

    #----Peak Intensity and Flux in 1cm Circle----
    peak_intensity=peak_density*AtomVelocity  #units of cm^-2s^-1

    flux = peak_intensity * flux_integral
    # sigma = 0
    # gamma = 2.58
    # spatial_gamma=2*gamma/10  #this is the FWHM in cm NOT HWHM in mm
    # spatial_sigma = sigma/10
    #
    # def integrand(z):
    #     return z*(spatial_gamma/2)**2/(z**2+(spatial_gamma/2)**2)
    #
    # radius=0.5 #cm
    # flux_integral = quad(integrand,0, radius)[0]
    # print('flux',2*np.pi * flux_integral)
    # flux=peak_intensity*2*np.pi*flux_integral

    if truth == True:

        print()
        print('|---Image and Box Location for Max Signal---|')
        print('|Max Pixel Value:', max_pixel_value)
        print('|Image on Fits File with Max Signal:', imageMax+1)
        print('|x pixel coordinates:', xMax, xMax+boxSize)
        print('|y pixel coordinates:', yMax, yMax+boxSize)
        print('|z1', z1Bound)
        print('|z2', z2Bound)
        print('|Max Total Counts:', countsMax)
        print()

        # print('Spatial gamma, FWHM (mm):',2*spatial_params[3])
        # print('Spatial FWHM (mm)',spatial_FWHM)

        if showPlots:

            # test=np.column_stack([spatial_y, spatial_signal])
            # data=np.column_stack([test,qG(spatial_y, *qGparams) ])
            # datafile = 'run30_spatialProfile.txt'
            # np.savetxt(datafile, data)
            print('spatial FWHM', spatialFWHM)

            plt.scatter(spatial_y, spatial_signal, c='r', label='data')
            plt.plot(spatial_y, qG(spatial_y, *qGparams), label='fit')
            plt.suptitle('Spatial Profile')
            #plt.title('Voigt FWHM ='+str(np.round(spatial_fwhm, 2))+'mm')
            plt.ylabel('Intensity (a.u.)',fontsize = 12)
            plt.xlabel('Transverse Position (mm)', fontsize = 12)
            plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
            #plt.grid()
            plt.legend()
            plt.show()

        #---Uncertainty in Spatial FWHM---

        spatial_fit = qG(spatial_y, *qGparams)
        # FWHM_fit_uncertainty = voigtFWHM_uncertainty(spatial_y, spatial_signal, spatial_fit, spatial_params, 0.001)
        FWHM_fit_uncertainty = qGaussFWHM_uncertainty(spatial_y, spatial_signal, spatial_fit, qGparams)
        print('FWHM uncertainty from fit',FWHM_fit_uncertainty)

        FWHM_uncertainty = np.sqrt((spatial_FWHM/magnification)**2*M_uncertainty**2 + FWHM_fit_uncertainty[0]**2)
        print('spatial FWHM', spatialFWHM)
        print('FWHM uncertainty including optics', FWHM_uncertainty)

        #---End of Spatial Uncertainty---

        #---Cross Section Uncertainty---
        #Correction factor for best estimate of uncertainty of voltage between the dips
        correctionFactor = 1/1.027
        area_under_curve=np.trapz(freqSignal, x=2*np.pi*imageFreqMhzArr*10**6*correctionFactor)
        maxSignal=np.max(freqSignal)

        S_w_test=maxSignal/area_under_curve  #s/rad
        S_w_uncertainty = np.abs(S_w-S_w_test)
        cross_section_test=(3/4)*gv.lambda_D2line**2*gv.A21*S_w_test
        cross_section_uncertainty = np.abs(cross_section-cross_section_test)

        #---End of cross section uncertainty---

        freqSignal = freqSignal / (np.max(freqSignal))

        #---Spectral Profiles---
        spectralFit=fit_Spectral_Data(imageFreqMhzArr, freqSignal, lensHeating=lensHeating, vTMaxLens=.043*211
                                      , peakMode='single', laserJitter=laserJitter)
        spectralFit.print_Results()

        spectralResults=fit_Spectral_Data(imageFreqMhzArr, freqSignal, lensHeating=False, vTMaxLens=.07*211*.9
                                      , peakMode='single', laserJitter=laserJitter)

        spectralResults.print_Results()

        spectralResults_LensHeating = fit_Spectral_Data_ConstrainedSigma(imageFreqMhzArr, freqSignal, 7, lensHeating=lensHeating, vTMaxLens=.043*211
                                      , peakMode='single', laserJitter=laserJitter)

        spectralResults_LensHeating.print_Results()

        if showPlots:
            plt.scatter(imageFreqMhzArr, freqSignal, c='r', label='Data')
            plt.plot(imageFreqMhzArr, spectralResults.fit_Result_Function(imageFreqMhzArr), label='Fit')
            plt.suptitle('Summed Signal vs. Image w/o Lens Heating')
            plt.xlabel('Relative Frequency (MHz)')
            plt.ylabel('Intensity (a.u.)')
            plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
            plt.legend()
            plt.show()

            plt.scatter(imageFreqMhzArr, freqSignal, c='r',label='Data')
            plt.plot(imageFreqMhzArr,spectralResults_LensHeating.fit_Result_Function(imageFreqMhzArr),label='Fit')
            plt.suptitle('Summed Signal vs. Image with Lens Heating')
            plt.xlabel('Relative Frequency (MHz)')
            plt.ylabel('Intensity (a.u.)')
            plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
            plt.legend()
            plt.show()


        #---Density Uncertainty---
        dz1 = z_uncertainty(magnification,binning,xMax,laser_profile_center,
                            pixel_size_on_camera,M_uncertainty,laser_center_uncertainty)


        dy = binning * boxSize * pixel_size_on_camera * M_uncertainty
        dwidth = dy
        y0 = pixelSize * boxSize

        dI = I_uncertainty(total_laser_power,y0,z1Bound,z2Bound,wz,power_uncertainty,dy,dz1,dwidth,wz_uncertainty)

        df = f_uncertainty(lens_diameter,lens_to_laser_length,diameter_uncertainty,laser_length_uncertainty)

        density_uncertainty = FluorDensityUncertainty(const,countsMax,f,eff,I,cross_section,
                                                      countsMax*0.01,df,eff_uncertainty,dI,cross_section_uncertainty) / 1e6
        #---End of Density Uncertainty---

        #---Intensity Uncertainty---
        peak_intensity_uncertainty=np.sqrt(
            AtomVelocity**2*density_uncertainty**2+peak_density**2*AtomVelocity_uncertainty**2)

        #---Flux Uncertainty---
        #Units need to be in cm

        sigma_uncertainty = FWHM_fit_uncertainty[1]
        q_uncertainty = FWHM_fit_uncertainty[2]
        eta_uncertainty = FWHM_fit_uncertainty[3]

        flux_integral_list = []

        print('test')
        for i in range(0,500):
            sigma_temp=np.random.normal(qGparams[3], sigma_uncertainty)
            q_temp=np.random.normal(qGparams[4], q_uncertainty)
            eta_temp=np.random.normal(qGparams[5], eta_uncertainty)

            left_FitValues=qG(left_x, 0, 1.0, 0, sigma_temp, q_temp, eta_temp)
            right_FitValues=qG(right_x, 0, 1.0, 0, sigma_temp, q_temp, eta_temp)

            left_integrand=np.abs(np.multiply(left_x, left_FitValues))
            right_integrand=np.abs(np.multiply(right_x, right_FitValues))
            left_flux_integral=np.trapz(left_integrand, left_x)
            right_flux_integral=np.trapz(right_integrand, right_x)

            flux_integral_temp=np.pi*(left_flux_integral+right_flux_integral)/100  #factor of 100 to convert mm^2 to cm^2
            flux_integral_list.append(flux_integral_temp)

        plt.hist(flux_integral_list,bins=20)
        plt.title('Histogram for Flux integral uncertainty')
        plt.show()
        flux_integral_uncertainty = np.std(flux_integral_list)

        fluxUncertainty = np.sqrt((flux_integral*peak_intensity_uncertainty)**2 +
                                  (peak_intensity * flux_integral_uncertainty)**2)
        # fluxUncertainty = delta_flux(peak_density,AtomVelocity,flux_integral,spatial_gamma/10,
        #                              density_uncertainty,FWHM_uncertainty/10,300)

        #---End of Flux Uncertainty---

        #---Brightness Calculation---

        guess=[np.argmax(freqSignal), np.max(freqSignal)-np.min(freqSignal), 9, 4, np.min(freqSignal)]
        spectralParams, pcov=spo.curve_fit(voigtFit, imageFreqMhzArr, freqSignal, p0=guess)

        spectralGamma=spectralParams[3]
        spectralSigma=spectralParams[2]
        freqSignal_fit=voigtFit(imageFreqMhzArr, *spectralParams)

        spectralFWHM=0.5346*2*spectralGamma+np.sqrt(0.2166*4*spectralGamma**2+2.35**2*spectralSigma**2)
        spectralFWHM_uncertainty=voigtFWHM_uncertainty(imageFreqMhzArr, freqSignal, freqSignal_fit, spectralParams,
                                                       0.01)

        V_transervse_spread=gv.cLight*spectralFWHM*10**6*2*np.pi/gv.w21

        V_transervse_spread_uncertainty=gv.cLight*spectralFWHM_uncertainty[0]*10**6*2*np.pi/gv.w21

        solid_angle=(V_transervse_spread/(AtomVelocity/100))**2

        peak_brightness=peak_intensity/solid_angle*10000  #factor of 10,000 to convert cm^-2 to m^-2

        peak_brightness_uncertainty = brightness_uncertainty(peak_intensity*10000, AtomVelocity/100, V_transervse_spread,
                                                             peak_intensity_uncertainty*10000, AtomVelocity_uncertainty/100,
                                                             V_transervse_spread_uncertainty)

        print()
        print('|-----Summary of Results-----|')
        print('|S(w0) (s/rad)', S_w)
        print('|S(w0) uncertainty (s/rad)', S_w_uncertainty)
        print('|Photons per second:', "{:e}".format(photons_per_second))
        print('|Laser Integration (Wm):', "{:e}".format(I))
        print('|Cross Section (m^2)', "{:e}".format(cross_section))
        print('|Cross Section Uncertainty (m^2):',"{:e}".format(cross_section_uncertainty))
        print()
        print('|Spectral FWHM (MHz):', spectralFWHM)
        print('|Spectral FWHM uncertainty (MHz):', spectralFWHM_uncertainty)
        print('|Transerve velosity spread FWHM (m/s):', V_transervse_spread)
        print('|Transverse velocity spread FWHM uncertainty (m/s):', V_transervse_spread_uncertainty)
        print()
        print('|Spatial FWHM (mm)):', spatial_FWHM)
        print('|Spatial FWHM uncertainty (mm):', FWHM_uncertainty)
        # print('|Spatial gamma and sigma (cm):', spatial_gamma, spatial_sigma)
        print()
        print('|Peak Density (cm^-3):', "{:e}".format(peak_density))
        print('|Density Uncertainty (cm^-3):',"{:e}".format(density_uncertainty))
        print()
        print('|Peak Intensity (cm^-2s^-1):',"{:e}".format(peak_intensity))
        print('|Peak Intensity uncertainty (cm^-2s^-1):', "{:e}".format(peak_intensity_uncertainty))
        print('|Flux in 1cm circle (atoms/s):', "{:e}".format(flux))
        print('|Flux uncertainty in 1cm circle (atoms/s):',"{:e}".format(fluxUncertainty))
        print('|Brightness (atoms s^-1 m^-2 sr^-1):', "{:e}".format(peak_brightness))
        print('|Brightness uncertainty (atoms s^-1 m^-2 sr^-1):', "{:e}".format(peak_brightness_uncertainty))

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