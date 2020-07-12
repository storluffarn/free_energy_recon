#*-encoding:utf-8  -*
from numpy import *
from itertools import cycle
import matplotlib.pyplot as plt
from sys import argv
from py_utils import *
K = 2

x = loadtxt(argv[1]+"/x_out.dat")
support = loadtxt(argv[1]+"/lambda_out.dat")


s,z = averge_slip_force(x,support,K,4*10**-10)

ax = plt.subplot(1,1,1)
for k in range(len(x[0,:])):
    ax.plot(support, -K*(x[:,k]-support), color="red", alpha=0.4)
ax.plot(s,z, "o", ms=10)
plt.show()
