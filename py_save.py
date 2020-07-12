#*-encoding:utf-8  -*
from numpy import *
from itertools import cycle
import matplotlib.pyplot as plt
from sys import argv
from py_utils import *
from scipy.optimize import curve_fit
kb = 1.380*10**-23
beta = (300*kb)**-1
sigma = 0.3*10**(-9)
K = 2

def f_fit(x, a,b,c):
    return a - b*abs(log(x/c))**(2.0/3.0)
def U(F_star):
    ## in units of kt
    return beta*(sigma*F_star/(2*pi))


data = loadtxt("f_av.txt")
v = data[:,0];f_av = data[:,1]
popt, pcov = curve_fit(f_fit, v, f_av)
U_mes = U(popt[0])

plt.plot(v, f_av)
plt.show()


