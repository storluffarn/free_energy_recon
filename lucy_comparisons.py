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
#mpl.rcParams["backend"] = "TkAgg"
import matplotlib.pyplot as plt
#font1 = {"family": "serif",
#         "size"  : "22",
#         "weight": "bold"}
#mpl.rc('text', usetex=True)
#mpl.rc("font", family="serif") 
import time
from local_lucy import local_lucy,lucy_max,lucy_summed_max, local_lucy2, local_lucy3, write_signal

E = 7
runs = 10
Z = 10
T = 1
Z_min = 2
Z_max = 8

system("./run -M {} -Z {} -E {} -N {} -d signals -t {}".format(runs,Z, E,2*T,T))
x_f = loadtxt("signals/x_out.dat")

traces = len(x_f[0,:])
samples = len(x_f[:,0])

LR_bins = int(((Z_max-Z_min)/10.0)*samples*0.5)
wham_bins = (Z_max-Z_min)*20

system("./run -M {} -Z {} -E {} -N {} -d signals -t {}".format(runs,Z, E,2*T,T))
x_r = loadtxt("signals/x_out.dat")
x = hstack((x_f, x_r))        # put them together for fair comparisons


support = loadtxt("signals/lambda_out.dat")

ITER = 25;
ITER_LOCAL = 20;
slack = 1.0
max_force = (2*pi*E*(beta**-1)/sigma)

signal_f, bin_centers = get_signal(x, support, LR_bins, Z_min, Z_max)
savetxt('signals/signal.dat',signal_f)
#signal_f = loadtxt("signals/sinedata.dat")
#signal_bi, bin_centers = get_signal_bi(x_f, x_r, support, LR_bins, Z_min, Z_max, traces, samples)
g_wham, z_valid = wham(x, wham_bins, K, len(x[0,:]), samples, support, Z_min, Z_max)
g_wham = set_constant(g_wham, z_valid, E)
#g_wham_bi, z_valid_bi = wham_bi(x_f, x_r, support, wham_bins, Z_min, Z_max, traces, samples, Z)
#g_wham_bi = set_constant(g_wham_bi, z_valid_bi, E)


fig = plt.figure()
#ax1 = plt.subplot(3,1,1)       # stop aggregate_max

#g_bi,z= lucy_max(signal_bi, 1, bin_centers,K, beta, ITER,  max_force*slack, (Z_min,Z_max))
#g_bi = set_constant(g_bi, z, E)

#g_f,z= lucy_max(signal_f, 1, bin_centers,K, beta, ITER,  max_force*slack, (Z_min,Z_max))
#g_f = set_constant(g_f, z, E)

#ax1.plot(z/sigma,beta*g_f, color="blue", label=r"$g^{LR}$")
#ax1.plot(z/sigma,beta*g_bi, color="red", label=r"$g^{LR}_{bi}$")
#ax1.plot(z_valid/sigma, beta*g_wham,"o", color="cyan", label=r"$g_{wham}$")
#ax1.plot(z_valid_bi/sigma, beta*g_wham_bi,"o", color="yellow", label=r"$g^{bi}_{wham}$")

#ax1.plot(z/sigma, beta*surface_potential(z, E), "--", color="black")
#ax1.grid()
#ax1.legend()
#ax1.set_xlabel(r"$z[\sigma]$",size=25)
#ax1.set_ylabel(r"$\Delta F_{surf} [k_{B}T]$",size=25)

#ax1.set_title(r"Simple max.")
ax2 = plt.subplot(3,1,2)       # stop aggregate_max

#g_bi,z= lucy_summed_max(signal_bi, 1, bin_centers,K, beta, ITER,  max_force*slack, (Z_min,Z_max))
#g_bi = set_constant(g_bi, z, E)

g_f,z= lucy_summed_max(signal_f, 1, bin_centers,K, beta, ITER,  max_force*slack, (Z_min,Z_max))

savetxt('signals/oraw.dat',g_f)

g_f = set_constant(g_f, z, E)

#write_signal(signal_f)

#for i, x in enumerate(z) :
#     g_f[i] = g_f[i] - (float(i)/1) / len(g_f)*x/(1*10**10)

ax2.plot(z/sigma,beta*g_f, color="blue", label=r"$g^{LR}$")
#ax2.plot(z/sigma,beta*g_bi, color="red", label=r"$g^{LR}_{bi}$")
ax2.plot(z_valid/sigma, beta*g_wham,"o", color="cyan", label=r"$g_{wham}$")
#ax2.plot(z_valid_bi/sigma, beta*g_wham_bi,"o", color="yellow", label=r"$g^{bi}_{wham}$")

ax2.plot(z/sigma, beta*surface_potential(z, E), "--", color="black")
ax2.grid()
#ax2.legend()
#ax2.set_xlabel(r"$z[\sigma]$",size=25)
ax2.set_ylabel(r"$\Delta F_{surf} [k_{B}T]$")

ax2.set_title(r"Summed max.")
#ax3 = plt.subplot(3,1,3)       # stop aggregate_max

#g_bi,z= local_lucy2(signal_bi, 1, bin_centers,K, beta, ITER_LOCAL,  max_force*slack)
#g_bi,z= local_lucy3(signal_bi, 1, bin_centers,K, beta, ITER_LOCAL, 2*E*(beta**-1))
#g_bi = set_constant(g_bi, z, E)

#g_f,z= local_lucy2(signal_f, 1, bin_centers,K, beta, ITER_LOCAL,  max_force*slack)
#g_f,z= local_lucy3(signal_f, 1, bin_centers,K, beta, ITER_LOCAL,  2*E*(beta**-1))
#g_f = set_constant(g_f, z, E)

#ax3.plot(z/sigma,beta*g_f, color="blue", label=r"$g^{LR}$")
#ax3.plot(z/sigma,beta*g_bi, color="red", label=r"$g^{LR}_{bi}$")
#ax3.plot(z_valid/sigma, beta*g_wham,"o", color="cyan", label=r"$g_{wham}$")
#ax3.plot(z_valid_bi/sigma, beta*g_wham_bi,"o", color="yellow", label=r"$g^{bi}_{wham}$")

#ax3.plot(z/sigma, beta*surface_potential(z, E), "--", color="black")
#ax3.grid()
#ax3.legend()
#ax3.set_xlabel(r"$z[\sigma]$",size=25)
#ax3.set_ylabel(r"$\Delta F_{surf} [k_{B}T]$",size=25)

#ax3.set_title(r"Local likelihood.")

fig.suptitle("E: {}, Z: {}, T: {}, Zmin: {}, Zmax: {}, runs: {}, \n ITER: {}, ITER-LOCAL: {}, slack {}".format(E,Z,T,Z_min,Z_max,runs,ITER,ITER_LOCAL,slack))

timestr = time.strftime("%Y%m%d-%H%M%S")

filename = 'plots/LR_summed_max_' + timestr

fig.savefig(filename)

