#*-encoding:utf-8  -*
from numpy import *
from itertools import cycle
from sys import argv
from py_utils import *
from scipy.optimize import curve_fit
from os import system
from py_deconv import *
import matplotlib as mpl
mpl.rcParams["backend"] = "TkAgg"
import matplotlib.pyplot as plt
font1 = {"family": "serif",
         "size"  : "22",
         "weight": "bold"}
mpl.rc('text', usetex=True)
mpl.rc("font", family="serif") 



# it might be the speed that is very different
E = 3
fig = plt.figure()
ax = plt.subplot(1,1,1)
system("./run -M 50 -Z 10 -N 5 -E {}  -d temp_data -t 0.5".format(E))
x = loadtxt("temp_data/x_out.dat")
support = loadtxt("temp_data/lambda_out.dat")

traces = len(x[0,:])
samples = len(x[:,0])
wham_bins=120
z_min=2;z_max=5;

g_wham, z_valid = wham(x, wham_bins, K, len(x[0,:]), samples, support, z_min, z_max)
g_wham = set_constant(g_wham, z_valid, E)

z = linspace(z_min,z_max, 1000)*sigma






ax = plt.subplot(1,1,1)
ax.plot(z_valid/sigma, beta*g_wham,"o", color="cyan", label=r"$g_{wham}$")
#ax.plot(z_valid_bi/sigma, beta*g_wham_bi,"o", color="yellow", label=r"$g^{bi}_{wham}$")
ax.plot(z/sigma, beta*surface_potential(z, E), "--", color="black")
ax.grid()
ax.legend()
ax.set_xlabel(r"$z[\sigma]$",size=25)
ax.set_ylabel(r"$\Delta F_{surf} [k_{B}T]$",size=25)

plt.show()








# axlist = []
# count=1
# for N in range(2,11):
#     axlist.append(plt.subplot(3,3,count))
#     signal, bin_centers = get_signal(x, support, N*300, 2, 8)
#     axlist[count-1].plot(bin_centers/sigma, -log(signal), label=r"$N = {}$".format(N*300))
#     axlist[count-1].grid()
#     axlist[count-1].set_xlabel(r"$z[\sigma]$",size=25)
#     axlist[count-1].set_ylabel(r"$\Delta F_{tot} [k_{B}T]$",size=25)
#     count+=1

# plt.show()

# # We are just going to look at some signals to see why it oscillates like it does
