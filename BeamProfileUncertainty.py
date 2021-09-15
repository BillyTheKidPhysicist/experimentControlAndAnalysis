import numpy as np
import scipy.special as sps
import scipy.optimize as spo
import matplotlib.pyplot as plt
from scipy import integrate
import os
from UncertaintyAnalysis_Functions import find_Chi_Squared_Uncertainty

#---Beam Waist Uncertainty---
path = "C:\Data\Runs\\farFIeldVerticalLaserImageCalibration"
os.chdir(path)
txtfile = open("7_19_2021_laserChar_t.txt", "r")

t_waistArr = np.loadtxt(txtfile)
txtfile.close()

xvalue_t = t_waistArr[:,0]
yvalue_t = t_waistArr[:,1]


def gaussian(x,x0,A,w):
    return A*np.exp(-2*(x-x0)**2/(w**2))

guess = (1.4,1.74,0.89)
params_t = spo.curve_fit(gaussian,xvalue_t,yvalue_t,p0=guess)[0]

print(params_t)
plt.scatter(xvalue_t,yvalue_t)
plt.plot(xvalue_t,gaussian(xvalue_t,*params_t))
plt.show()

optimal_w=params_t[2]
print('Best estimate for w_t:',optimal_w)
test_size=0.075
chi_squared=[]
w_test=[]
for w_value in np.arange(optimal_w-test_size, optimal_w+test_size, test_size/50):

    w_test.append(w_value)


    def gaussianFixed(x, x0, A):
        w=w_value
        return gaussian(x,x0,A,w)


    guess_fixed=[1.4,1.7]

    params_t_fixed=spo.curve_fit(gaussianFixed, xvalue_t, yvalue_t,
                                       p0=guess_fixed)[0]

    R2=[]
    for i in np.arange(0, len(xvalue_t)):
        value_fit=gaussianFixed(xvalue_t[i], *params_t_fixed)
        value_data=yvalue_t[i]
        difference=(value_fit-value_data)**2
        R2.append(difference)

    chi_squared.append(np.sum(R2))

    # plt.plot(xvalue_t,gaussianFixed(xvalue_t,*params_t_fixed))
    # plt.scatter(xvalue_t,yvalue_t)
    # plt.show()

find_Chi_Squared_Uncertainty(w_test, chi_squared)

txtfile = open("7_20_2021_laserChar_l.txt", "r")

l_waistArr = np.loadtxt(txtfile)
txtfile.close()

xvalue_l = l_waistArr[:,0]
yvalue_l = l_waistArr[:,1]


guess = (60,10,28)
params_l = spo.curve_fit(gaussian,xvalue_l,yvalue_l,p0=guess)[0]

print(params_l)
plt.scatter(xvalue_l,yvalue_l)
plt.plot(xvalue_l,gaussian(xvalue_l,*params_l))
plt.show()

optimal_w=params_l[2]
print('Best estimate for w_l:',optimal_w)
test_size=4.5
chi_squared=[]
w_test=[]
for w_value in np.arange(optimal_w-test_size, optimal_w+test_size, test_size/50):

    w_test.append(w_value)


    def gaussianFixed(x, x0, A):
        w=w_value
        return gaussian(x,x0,A,w)


    guess_fixed=[60,10]

    params_t_fixed=spo.curve_fit(gaussianFixed, xvalue_l, yvalue_l,
                                       p0=guess_fixed)[0]

    R2=[]
    for i in np.arange(0, len(xvalue_l)):
        value_fit=gaussianFixed(xvalue_l[i], *params_t_fixed)
        value_data=yvalue_l[i]
        difference=(value_fit-value_data)**2
        R2.append(difference)

    chi_squared.append(np.sum(R2))

    # plt.plot(xvalue_t,gaussianFixed(xvalue_t,*params_t_fixed))
    # plt.scatter(xvalue_t,yvalue_t)
    # plt.show()

find_Chi_Squared_Uncertainty(w_test, chi_squared)
