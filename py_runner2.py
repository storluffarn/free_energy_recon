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




E = 7
system("./run -M 50 -Z 10 -E {}  -d temp_data".format(E))
x = loadtxt("temp_data/x_out.dat")
support = loadtxt("temp_data/lambda_out.dat")

f_av = average_slip_force2(x,support,K,4*10**-10)


signal, bin_centers = get_signal(x, support, 1000, 2, 8)
ITER = 200

g,z = LR_deconv(signal, 20, bin_centers,K, beta, ITER)
#g_stop,z = LR_deconv_maxforce(signal, 20, bin_centers,K, beta, ITER, f_av*2)
g_stop,z = LR_deconv_maxforce(signal, 20, bin_centers,K, beta, ITER, f_av)
g_stop_proper,z = LR_deconv_maxforce(signal, 20, bin_centers,K, beta, ITER, 2*pi*E*(beta**-1)/sigma )
g_stop_check,z =  LR_deconv_maxforce(signal, 20, bin_centers,K, beta, ITER, f_av*2)

g_stop= set_constant(g_stop, z, E)
g = set_constant(g, z, E)
g_stop_proper = set_constant(g_stop_proper, z, E)
g_stop_check = set_constant(g_stop_check, z, E)


ax = plt.subplot(2,1,1)
ax.plot(z/sigma,beta*g_stop, color="red", label=r"$g_{stop}$")
ax.plot(z/sigma,beta*g, color="blue", label=r"$g$")
ax.plot(z/sigma,beta*g_stop_proper, color="green", label=r"$g_{truestop}$")
ax.plot(z/sigma,beta*g_stop_check, color="cyan", label=r"$g_{check}$")

ax.plot(z/sigma, beta*surface_potential(z, E), "--", color="black")
ax.grid()
ax.legend()
ax.set_xlabel(r"$z[\sigma]$",size=25)
ax.set_ylabel(r"$\Delta F_{surf} [k_{B}T]$",size=25)
ax.set_title(r"$U_{{0}} = {} [k_{{B}}T]$".format(E), size=30)
#ax.set_ylabel(r"$\Delta F_{tot} [k_{B}T]$")


ax2 = plt.subplot(2,1,2)
ax2.plot(bin_centers/sigma, -log(signal))
ax2.grid()
ax2.set_xlabel(r"$z[\sigma]$",size=25)
ax2.set_ylabel(r"$\Delta F_{tot} [k_{B}T]$",size=25)


plt.show()
#signal, bin_centers, bin_edges, weight_work,support,ext_work= get_signal(traces, n_scans, scan_length, d_lambda,bins, recon_span)
