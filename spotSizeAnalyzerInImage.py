import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import scipy.special as sps

import scipy.optimize as spo

fileName='run15ShutterClosed'
fitsFile = fits.open(fileName+'.fits')
imagesArr = fitsFile[0].data
imagesArr=imagesArr.astype(float)
imagesArr=np.flip(imagesArr,axis=1)
trimUp=1400 #!!!!THIS IS SENSITIVE
imagesArr[imagesArr>trimUp]=trimUp
trimLow=1000
imagesArr[imagesArr<trimLow]=trimLow
imageSum=np.mean(imagesArr[20:30],axis=0)[150:-85,260:-230]
#print(imageSum.shape)
# plt.imshow(imageSum)
# plt.show()
y=np.mean(imageSum,axis=1)
x=np.linspace(0,0.2*y.shape[0],num=y.shape[0])

def fit(x,x0,sigma,gamma,a,b):
    z0=sps.voigt_profile(0,sigma,gamma)
    z=sps.voigt_profile(x-x0,sigma,gamma)
    return a*z/z0 +b
print(imagesArr.shape)
test=imagesArr[:,150:-85,260:-230]
print(test.shape)
test=np.mean(test,axis=2)[:,125]
print(test)
plt.plot(test)
plt.show()






guess=[25,1,1,100,1150]
bounds=(0,np.inf)
params,pcov=spo.curve_fit(fit,x,y,p0=guess,bounds=bounds)
perr = np.sqrt(np.diag(pcov))
fL=params[2]*2
fG=2*params[1]*np.sqrt(2*np.log(2))
fV=.5346*fL+np.sqrt(.2166*fL**2+fG**2)
print(fV,fG,fL)
xTest=np.linspace(x[0],x[-1],num=1000)
yTest=fit(xTest,*params)-params[-1]


plt.scatter(x,y-params[-1],marker='x',label='Data')
plt.plot(xTest,yTest,c='r',label='Fit')
plt.legend()
plt.title('Transverse signal fitted with voigt')
plt.grid()
plt.ylabel('Arbitrary value')
plt.xlabel('Millimeters')
plt.savefig('plot')
plt.show()
