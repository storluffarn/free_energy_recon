#*-encoding:utf-8  -*
from numpy import *
from itertools import cycle
#import matplotlib.pyplot as plt
from sys import argv
from py_utils import *
from py_deconv import *
from os import system
from scipy.optimize import curve_fit
import matplotlib as mpl
mpl.rcParams["backend"] = "TkAgg"
import matplotlib.pyplot as plt
font1 = {"family": "serif",
         "size"  : "22",
         "weight": "bold"}
mpl.rc('text', usetex=True)
mpl.rc("font", family="serif") 

# so it scales perfectly

E = 13
system("./run -M 10 -Z 2 -E {} -t {} -N {} -d double".format(E,1,10))
x = loadtxt("double/x_out.dat")
support = loadtxt("double/lambda_out.dat")
traces = len(x[0,:])
samples = len(x[:,0])

wham_bins=40
Z_min = 0;Z_max=1.5

signal, bin_centers = get_signal(x, support, 300, Z_min, Z_max)
ITER = 1000
g,z = LR_deconv(signal, 20, bin_centers,K, beta, ITER)
g = set_constant(g, z, E)

g_wham, z_valid = wham(x, wham_bins, K, len(x[0,:]), samples, support, Z_min, Z_max)
g_wham = set_constant(g_wham, z_valid, E)


ax = plt.subplot(2,1,1)
ax.plot(z/sigma,beta*g, color="blue", label=r"$g$")
ax.plot(z/sigma, beta*surface_potential(z, E), "--", color="black")
ax.plot(z_valid/sigma, beta*g_wham, "o", label=r"g_{wham}")
ax.grid()
ax.legend()
ax.set_xlim([-0.5, 1.5])
ax.set_ylim([-0.3*E, E*2])
ax.set_xlabel(r"$z[\sigma]$",size=25)
ax.set_ylabel(r"$\Delta F_{surf} [k_{B}T]$",size=25)
ax.set_title(r"$U_{{0}} = {} [k_{{B}}T]$".format(E), size=30)

ax2 = plt.subplot(2,1,2)
ax2.plot(bin_centers/sigma, -log(signal))
ax2.grid()
ax2.set_xlabel(r"$z[\sigma]$",size=25)
ax2.set_ylabel(r"$\Delta F_{tot} [k_{B}T]$",size=25)


plt.show()



# the problem with this guy seems to be that it UNDERestimates, no matter how many iterations you do
# for k in range(10):
#     plt.plot(support, x[:,k])
