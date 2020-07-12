#*-encoding:utf-8  -*
from numpy import *
from itertools import cycle
import matplotlib.pyplot as plt
from sys import argv
from py_utils import *
from scipy.optimize import curve_fit
import numpy.polynomial.polynomial as pol
kb = 1.380*10**-23
beta = (300*kb)**-1
sigma = 0.3*10**(-9)
K = 2

E_0 = 64*(beta**-1)*(sigma**-4)
E = E_0*10


# lets just try to plot it and see what we get...
# might be good


def pol(x):
    return E*(  0.25*x**(4) - 0.5*sigma*x**(3) + 0.25*(sigma**2)*(x**2)   )



t = linspace(-0.5*sigma, 1.5*sigma, 1000)
y = pol(t)

plt.plot(t/sigma,y*beta)
plt.grid()
plt.show()
