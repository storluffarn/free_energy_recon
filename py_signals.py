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

from local_lucy import local_lucy

# task... check if the bidirectional estimator produces better result. It would e.g. be terrific if it removed the bias in the data
E=7;
# traces = 10
# samples = 5000

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
#wham_bins = int(80*(float(samples)/500))
wham_bins = (Z_max-Z_min)*20

system("./run -M {} -Z {} -E {} -N {} -d signals -t {}".format(runs,Z, E,2*T,T))
x_r = loadtxt("signals/x_out.dat")

x = hstack((x_f, x_r))        # put them together for fair comparisons
 
support = loadtxt("signals/lambda_out.dat")

signal_f, bin_centers = get_signal(x, support, LR_bins, Z_min, Z_max)
signal_bi, bin_centers = get_signal_bi(x_f, x_r, support, LR_bins, Z_min, Z_max, traces, samples)


ITER = 100
ITER_LOCAL=25
slack = 1.00
max_force = (2*pi*E*(beta**-1)/sigma)
#g_bi,z = LR_deconv(signal_bi, 20, bin_centers,K, beta, ITER)
g_bi,z= LR_deconv_maxforce(signal_bi, 1, bin_centers,K, beta, ITER,  max_force*slack, (Z_min,Z_max))
g_bi_local,z= local_lucy(signal_bi, 1, bin_centers,K, beta, ITER_LOCAL,  max_force*slack)

g_bi = set_constant(g_bi, z, E)
g_bi_local = set_constant(g_bi_local, z, E)

#g_f,z = LR_deconv(signal_f, 20, bin_centers,K, beta, ITER)
g_f,z= LR_deconv_maxforce(signal_f, 1, bin_centers,K, beta, ITER, max_force*slack, (Z_min, Z_max) )
g_f_local,z= local_lucy(signal_f, 1, bin_centers,K, beta, ITER_LOCAL, max_force*slack)
#g_f,z= LR_deconv_maxforce_pbc(signal_f, int(LR_bins*0.5/3.0), bin_centers,K, beta, ITER, 2*pi*E*(beta**-1)/sigma )
g_f = set_constant(g_f, z, E)
g_f_local = set_constant(g_f_local, z,E)


print "deconv window: ", int(LR_bins*0.5/3.0)

g_wham, z_valid = wham(x, wham_bins, K, len(x[0,:]), samples, support, Z_min, Z_max)
g_wham = set_constant(g_wham, z_valid, E)

g_wham_bi, z_valid_bi = wham_bi(x_f, x_r, support, wham_bins, Z_min, Z_max, traces, samples, Z)
g_wham_bi = set_constant(g_wham_bi, z_valid_bi, E)


ax = plt.subplot(2,1,1)
ax.set_title(r"$(U_{{0}}, s, N,N^{{LR}}_{{bins}},N^{{wham}}_{{bins}} ) = ({} [k_{{B}}T], {},{},{},{})$".format(E, samples, runs*2, LR_bins, wham_bins), size=30)
ax.plot(z/sigma,beta*g_f, color="blue", label=r"$g^{LR}$")
ax.plot(z/sigma,beta*g_bi, color="red", label=r"$g^{LR}_{bi}$")

ax.plot(z/sigma,beta*g_f_local,"--", color="blue", label=r"$g^{LR}(local)$")
ax.plot(z/sigma,beta*g_bi_local,"--" ,color="red", label=r"$g^{LR}_{bi}(local)$")


ax.plot(z_valid/sigma, beta*g_wham,"o", color="cyan", label=r"$g_{wham}$")
ax.plot(z_valid_bi/sigma, beta*g_wham_bi,"o", color="yellow", label=r"$g^{bi}_{wham}$")

ax.plot(z/sigma, beta*surface_potential(z, E), "--", color="black")
ax.grid()
ax.legend()
ax.set_xlabel(r"$z[\sigma]$",size=25)
ax.set_ylabel(r"$\Delta F_{surf} [k_{B}T]$",size=25)

ax2 = plt.subplot(2,1,2)
ax2.plot(bin_centers/sigma, -log(signal_f), label=r"Forward dynamics.")
ax2.plot(bin_centers/sigma, -log(signal_bi), label=r"Bidirectional.")
ax2.grid()
ax2.legend()
ax2.set_xlabel(r"$z[\sigma]$",size=25)
ax2.set_ylabel(r"$\Delta F_{tot} [k_{B}T]$",size=25)



plt.show()
