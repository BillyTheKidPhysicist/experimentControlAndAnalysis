import numpy as np
from astropy.io import fits
import os
import matplotlib.pyplot as plt
import scipy.optimize as spo
from scipy.integrate import quad
from UncertaintyAnalysis_Functions import find_Chi_Squared_Uncertainty
from UncertaintyAnalysis_Functions import delta_n
from UncertaintyAnalysis_Functions import delta_flux

# plt.rcParams["figure.figsize"]=(8, 6)

def load_Array_Of_Images_From_Fits(path,fileName):
    #this will open a fits file found at path+fileName and return an array of images. If the fits file contains a single
    #image it will still be return as an array of images, just one image. The return array has shape (n,x1,x2) where n
    #is the number of images, and x1 and x2 are image dimensions. Orientation has not been checked here.
    print('correct orientation is not guranteed. Double check')
    os.chdir(path)

    #Opening fits file and creating array. Flip the y axis such that image is oriented correctly.
    fitsFile=fits.open(fileName+'.fits')
    imageDataDimensionality=len(fitsFile[0].data.shape) #may be 2d or 3d!

    if imageDataDimensionality==3: #list of images
        imagesList=fitsFile[0].data
    elif imageDataDimensionality==2:
        imagesList=[fitsFile[0].data]
    else:
        raise Exception("image format is wrong")
    imagesArr=np.asarray(imagesList).astype(float) #convert from ints to floats!
    return imagesArr

def fit_Result_Function(x,x0,a,b,gamma):
    #able transformed absorption lorentzian
    #a: positive value, dip is this depth
    #gamma: FWHM
    return -a*(gamma/2)/np.sqrt((x-x0)**2+(gamma/2)**2) +b

path='C:\Data\Runs\9_23_21\\run15Folder'
fileName='run15_absorptionSignalImages'
imagesArr=load_Array_Of_Images_From_Fits(path,fileName)

# print(imagesArr.shape)

# imageList=[imageList[8]]

# for imageRaw in imagesArr[0:-1]:
for imageRaw in imagesArr[3:4]:
    binning=1
    magnification=3.05
    M_uncertainty = 0.02
    pixelSizeOnCamera=24e-6 #m
    pixelSizeEffective=pixelSizeOnCamera*binning*magnification
    yIndexStart=280
    yIndexEnd=620
    xIndexStart=550#600
    xIndexEnd=650#650

    #box dimensions
    y1 = 190
    y2 = 840
    x1 = 170
    x2 = 820
    image=imageRaw.copy() #to no modify original
    image2=np.flip(image[x1:x2,y1:y2],axis=0)
    p = np.linspace(0,image.shape[0]*pixelSizeEffective * 100,num=image.shape[0])
    #p = np.linspace(0, image2.shape[0], num=image2.shape[0])

    plt.imshow(image2,vmin=0.849,vmax=1.0,extent=[p[x1]-1.25,p[x2]-1.25,p[y1]-1.25,p[y2]-1.25])
    # plt.axhline(y=yIndexStart*pixelSizeEffective*100, c='r', linestyle="--")
    # plt.axhline(y=yIndexEnd*pixelSizeEffective*100, c='r', linestyle="--")
    # plt.axvline(x=xIndexStart*pixelSizeEffective*100, c='r', linestyle="--")
    # plt.axvline(x=xIndexEnd*pixelSizeEffective*100, c='r', linestyle="--")
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=12)
    cb.set_label(label='Transmittance',fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel('Transverse Position (cm)',fontsize=12)
    plt.xticks(fontsize=12)
    plt.xlabel('Longitudinal Position (cm)',fontsize=12)
    plt.savefig('absorptionImage.png', dpi = 1200)
    plt.show()
    assert xIndexStart<image.shape[1] and xIndexEnd<image.shape[1]
    assert yIndexStart<image.shape[0] and yIndexEnd<image.shape[0]

    yArr0=np.arange(image.shape[0])*pixelSizeEffective
    profile0=np.mean(image[:,xIndexStart:xIndexEnd],axis=1)

    yArr=yArr0[yIndexStart:yIndexEnd].copy()
    profile=profile0[yIndexStart:yIndexEnd].copy()

    yValueAtDip=yArr[np.argmin(profile)]
    yArr-=yValueAtDip #center near zero
    yArrplot=yArr + (yIndexEnd+yIndexStart)*pixelSizeEffective / 2

    #x0,a,b,gamma
    guess=[2.0e-5,.12,1.0,2.5e-3]
    params=spo.curve_fit(fit_Result_Function,yArr,profile,p0=guess)[0]
    print('Fit Parameters',params)
    FWHM = params[3]*1000
    plt.plot(yArr*10,profile,label='Data')
    plt.plot(yArr*10,fit_Result_Function(yArr,*params),label='Fit')
    plt.title('Transmission vs. Transverse Position')
    plt.ylabel('Percent Transmission')
    plt.xlabel('Transverse Position (cm)')
    plt.legend()
    plt.grid()
    plt.show()

    #-----Uncertainty in gamma and % absorption------

    optimal_a = params[1]
    chi_squared_a = []
    test_a = []
    idk = 0.85

    sigma_squared = np.square(np.subtract(profile,fit_Result_Function(yArr,*params)))

    for a_values in np.linspace(optimal_a*0.995,optimal_a*1.0015, num = 500 ):

        test_a.append(a_values)
        def fit_Result_fixed_a(x, x0, b, gamma):
            a = a_values
            return fit_Result_Function(x,x0,a,b,gamma)

        guess_a_fixed = [0,0.1,0.001]

        params_a_fixed = spo.curve_fit(fit_Result_fixed_a,yArr,profile,p0=guess_a_fixed)[0]
        # print(params_a_fixed)
        # print(a_values)
        # plt.plot(yArr,profile)
        #plt.title('fit results for test values of a')
        # plt.plot(yArr,fit_Result_fixed_a(yArr,*params_a_fixed),c='r')
        # plt.show()

        R_a=np.subtract(profile, fit_Result_fixed_a(yArr, *params_a_fixed))
        R2_a=np.square(R_a)

        chi_squared=np.sum(np.divide(R2_a, sigma_squared))

        chi_squared_a.append(chi_squared)

    print('--absorption--')
    print('Optimal a:',optimal_a)
    plt.scatter(test_a, chi_squared_a)
    plt.show()
    bounds_a = find_Chi_Squared_Uncertainty(test_a,chi_squared_a,True)
    a_uncertainty = bounds_a
    #a_uncertainty = max(np.abs(optimal_a-bounds_a[0]),np.abs(optimal_a-bounds_a[1]))
    print('a uncertainty:',a_uncertainty)


    optimal_gamma = params[3]
    chi_squared_gamma = []
    test_gamma = []
    idk = 0.55

    for gamma_values in np.linspace(optimal_gamma*0.995, optimal_gamma*1.005, num = 500):

        test_gamma.append(gamma_values)
        def fit_Result_fixed_gamma(x, x0, a, b):
            gamma = gamma_values
            return fit_Result_Function(x,x0,a,b,gamma)

        guess_gamma_fixed = [0,0.1,0.005]

        params_gamma_fixed = spo.curve_fit(fit_Result_fixed_gamma,yArr,profile,p0=guess_gamma_fixed)[0]

        # print(params_gamma_fixed)
        # print(gamma_values)
        # plt.plot(yArr,profile)
        # plt.plot(yArr,fit_Result_fixed_gamma(yArr,*params_gamma_fixed),c='r')
        # plt.show()

        R_gamma=np.subtract(profile, fit_Result_fixed_gamma(yArr, *params_gamma_fixed))
        R2_gamma=np.square(R_gamma)

        chi_squared=np.sum(np.divide(R2_gamma, sigma_squared))

        chi_squared_gamma.append(chi_squared)

    print('--gamma--')
    print('Optimal gamma:',optimal_gamma)
    bounds_gamma = find_Chi_Squared_Uncertainty(test_gamma,chi_squared_gamma,True)
    gamma_uncertainty = bounds_gamma
    #gamma_uncertainty=max(np.abs(optimal_gamma-bounds_gamma[0]), np.abs(optimal_gamma-bounds_gamma[1]))
    #Including uncertainty in magnification
    gamma_uncertainty = np.sqrt(gamma_uncertainty**2+(optimal_gamma*M_uncertainty/magnification)**2)

    print('gamma uncertainty [cm]:', gamma_uncertainty*100)

    #-----End of Uncertainty-----

    #-----Jeremy's addition to compute peak density-----

    #User needs to impute S(w) in s/rad to be able to compute the frequency integrated cross section.
    #User needs to input bounds of spatial integration. Will treat this as diameter of Obs Cell, so ~10cm
    #Don't expect S(w) to vary much so it will probably be around 6.4*10^-9 s/rad
    #Will treat spatial profile as Lorentzian and will use gamma (FWHM) of spatial fit of absorption data.

    S_w = 5.704886*10**-9
    cross_section_uncertainty = 1.917501e-15

    def cross_section(S):
        crossSection_in_m2 = 0.75 * (670.776*10**-9)**2*3.689*10**7*S * 2/3
        crossSection_in_cm2 = crossSection_in_m2 * 100**2
        return crossSection_in_cm2


    #ALL gamma should be FWHM NOT HWHM!!
    gamma = params[3]*100 #in units of cm
    integrationBound = 5 #in units of cm

    def spatial_profile_integrand (x):
        return (gamma/2)**2 / (x**2+(gamma/2)**2)

    integrated_spatial_profile = quad(spatial_profile_integrand, -integrationBound, integrationBound)

    absorption_ratio = fit_Result_Function(params[0],*params)
    print()
    print('Percent Absorbed', (1 - absorption_ratio)*100)

    peak_density = -np.log(absorption_ratio) / (cross_section(S_w) * integrated_spatial_profile[0])
    print()


    #input needs to be in cm!
    densityUncertainty=delta_n(absorption_ratio, cross_section(S_w), 5, gamma, integrated_spatial_profile[0],
                               gamma_uncertainty*100, cross_section_uncertainty*100*100,a_uncertainty)

    atomVelocity = 21100 #cm/s
    peak_intensity = peak_density * atomVelocity

    #---Computing Flux---

    def integrand(z):
        return z*(gamma/2)**2/(z**2+(gamma/2)**2)

    radius=0.5  #cm
    flux_integral=quad(integrand, 0, radius)[0]
    flux=peak_intensity*2*np.pi*flux_integral

    #Units of input need to be in cm
    fluxUncertainty = delta_flux(peak_density,atomVelocity,flux_integral,gamma,densityUncertainty,gamma_uncertainty,300)


    #---Optical Density---
    print('also a',absorption_ratio)
    print('a',optimal_a)
    OD = -np.log10(absorption_ratio)

    delta_OD = 1/((absorption_ratio)*np.log10(10))*a_uncertainty

    print()
    print('Optical Density:',OD)
    print('Optical Density Uncertainty:',delta_OD)
    print('FWHM (mm):', params[3]*1000)
    print('FWHM Uncertainty (mm):', gamma_uncertainty*1000)
    print('Peak Density (cm^-3):', "{:e}".format(peak_density))
    print('Peak Density Uncertainty (cm^-3):', "{:e}".format(densityUncertainty))
    print('Peak Intensity:',"{:e}".format(peak_intensity))
    print('Flux through 1cm circle:',"{:e}".format(flux))
    print('Flux Uncertainty:',"{:e}".format(fluxUncertainty))
    print()
    print('-----Next Image-----')
    print()

    # test=np.column_stack([yArrplot*1000, profile])
    # data=np.column_stack([test,fit_Result_Function(yArr, *params)])
    # datafile = 'run26_absorption_spatialProfile.txt'
    # np.savetxt(datafile, data)

    plt.plot(yArrplot*1000, profile, label='Data')
    plt.plot(yArrplot*1000, fit_Result_Function(yArr, *params), label='Fit')
    #plt.title('Transmission vs. Transverse Position')
    plt.ylabel('Transmission %',fontsize=20)
    plt.xlabel('Transverse Position (mm)',fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend()

    plt.show()