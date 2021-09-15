import numpy as np
from astropy.io import fits
import os
import matplotlib.pyplot as plt
import scipy.optimize as spo
from scipy.integrate import quad
from UncertaintyAnalysis_Functions import find_Chi_Squared_Uncertainty
from UncertaintyAnalysis_Functions import delta_n

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
    #able transformed absorption lorentzia
    #a: positive value, dip is this depth
    #gamma: FWHM
    return -a*(gamma/2)/np.sqrt((x-x0)**2+(gamma/2)**2) +b

path='C:\Data\Runs\9_6_21\\run12Folder'
fileName='run12_absorptionSignalImages'
imagesArr=load_Array_Of_Images_From_Fits(path,fileName)
print(imagesArr.shape)

# imageList=[imageList[8]]

for imageRaw in imagesArr:
    binning=1
    magnification=3.0
    pixelSizeOnCamera=24e-6 #m
    pixelSizeEffective=pixelSizeOnCamera*magnification*binning
    yIndexStart=300
    yIndexEnd=600
    xIndexStart=550
    xIndexEnd=650
    image=imageRaw.copy() #to no modify original
    poop = np.linspace(0,image.shape[0]*pixelSizeEffective * 100,num=image.shape[0])

    # image[image<.9]=.9
    # image[image>1.05]=1.05
    plt.imshow(image,vmin=0.88,vmax=1.05,extent=[poop[0],poop[-1],poop[0],poop[-1]])
    # plt.axhline(y=yIndexStart, c='r', linestyle="--")
    # plt.axhline(y=yIndexEnd, c='r', linestyle="--")
    # plt.axvline(x=xIndexStart, c='r', linestyle="--")
    # plt.axvline(x=xIndexEnd, c='r', linestyle="--")
    plt.colorbar()
    plt.ylabel('Transverse Position (cm)')
    plt.xlabel('Longitudinal Position (cm)')
    plt.show()
    assert xIndexStart<image.shape[1] and xIndexEnd<image.shape[1]
    assert yIndexStart<image.shape[0] and yIndexEnd<image.shape[0]

    yArr0=np.arange(image.shape[0])*pixelSizeEffective
    profile0=np.mean(image[:,xIndexStart:xIndexEnd],axis=1)

    yArr=yArr0[yIndexStart:yIndexEnd].copy()
    profile=profile0[yIndexStart:yIndexEnd].copy()

    yValueAtDip=yArr[np.argmin(profile)]
    yArr-=yValueAtDip #center near zero


    #x0,a,b,gamma
    guess=[0.0,.1,1.0,2.5e-3]
    params=spo.curve_fit(fit_Result_Function,yArr,profile,p0=guess)[0]
    print('Fit Parameters',params)
    print('FWHM [mm]',params[3]*1000)
    # plt.plot(yArr*10,profile,label='Data')
    # plt.plot(yArr*10,fit_Result_Function(yArr,*params),label='Fit')
    # plt.title('Transmission vs. Transverse Position')
    # plt.ylabel('Percent Transmission')
    # plt.xlabel('Transverse Position (cm)')
    # plt.legend()
    # plt.grid()
    # plt.show()

    #-----Uncertainty in gamma and % absorption------

    optimal_a = params[1]
    chi_squared_a = []
    test_a = []
    idk = 0.90

    for a_values in np.arange(optimal_a*idk,optimal_a*(2-idk),optimal_a*idk / 500 ):

        test_a.append(a_values)
        def fit_Result_fixed_a(x, x0, b, gamma):
            a = a_values
            return fit_Result_Function(x,x0,a,b,gamma)

        guess_a_fixed = [0,0.1,0.001]
        params_a_fixed = spo.curve_fit(fit_Result_fixed_a,yArr,profile,p0=guess_a_fixed)[0]
        # print(params_a_fixed)
        # print(a_values)
        # plt.plot(yArr,profile)
        # plt.plot(yArr,fit_Result_fixed_a(yArr,*params_a_fixed),c='r')
        # plt.show()
        R2_a=[]
        for i in np.arange(0, len(yArr)):
            value_fit=fit_Result_fixed_a(yArr[i], *params_a_fixed)
            value_data=profile[i]
            difference=(value_fit-value_data)**2
            R2_a.append(difference)

        chi_squared_a.append(np.sum(R2_a))

    print('absorption')
    find_Chi_Squared_Uncertainty(test_a,chi_squared_a)

    optimal_gamma = params[3]
    chi_squared_gamma = []
    test_gamma = []
    idk = 0.8
    for gamma_values in np.arange(optimal_gamma*idk, optimal_gamma*(2-0.9*idk), optimal_gamma*idk/200):

        test_gamma.append(gamma_values*100)
        def fit_Result_fixed_gamma(x, x0, a, b):
            gamma = gamma_values
            return fit_Result_Function(x,x0,a,b,gamma)

        guess_gamma_fixed = [0,0.01,0.005]
        params_gamma_fixed = spo.curve_fit(fit_Result_fixed_gamma,yArr,profile,p0=guess_gamma_fixed)[0]

        # print(params_gamma_fixed)
        # print(gamma_values)
        # plt.plot(yArr,profile)
        # plt.plot(yArr,fit_Result_fixed_gamma(yArr,*params_gamma_fixed),c='r')
        # plt.show()

        R2_gamma = []

        for i in np.arange(0, len(yArr)):
            value_fit=fit_Result_fixed_gamma(yArr[i], *params_gamma_fixed)
            value_data=profile[i]
            difference=(value_fit-value_data)**2
            R2_gamma.append(difference)

        chi_squared_gamma.append(np.sum(R2_gamma))

    print('gamma')
    find_Chi_Squared_Uncertainty(test_gamma,chi_squared_gamma)

    #-----End of Uncertainty-----

    #-----Jeremy's addition to compute peak density-----

    #User needs to impute S(w) in s/rad to be able to compute the frequency integrated cross section.
    #User needs to input bounds of spatial integration. Will treat this as diameter of Obs Cell, so ~10cm
    #Don't expect S(w) to vary much so it will probably be around 6.4*10^-9 s/rad
    #Will treat spatial profile as Lorentzian and will use gamma (FWHM) of spatial fit of absorption data.

    S_w = 6.03*10**-9

    def cross_section(S):
        crossSection_in_m2 = 0.75 * (670.776*10**-9)**2*3.689*10**7*S / 2
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
    test=delta_n(absorption_ratio, cross_section(S_w), 5, gamma, integrated_spatial_profile[0], 0.07, cross_section(S_w)*0.05,
                 0.005)

    print('Peak Density [cm^-3]:', "{:e}".format(peak_density))
    print('Peak Density Uncertainty [cm^-3]:',"{:e}".format(test))
    print()
    atomVelocity = 20600
    peak_intensity = peak_density * atomVelocity
    def integrand(z):
        return z*(gamma/2)**2/(z**2+(gamma/2)**2)


    radius=0.5  #cm
    flux=peak_intensity*2*np.pi*quad(integrand, 0, radius)[0]

    print('Peak Intensity:',"{:e}".format(peak_intensity))
    print('Flux through 1cm circle:',"{:e}".format(flux))

    plt.plot(yArr*10, profile, label='Data')
    plt.plot(yArr*10, fit_Result_Function(yArr, *params), label='Fit')
    plt.title('Transmission vs. Transverse Position')
    plt.ylabel('Percent Transmission')
    plt.xlabel('Transverse Position (cm)')
    plt.legend()
    plt.grid()
    plt.show()