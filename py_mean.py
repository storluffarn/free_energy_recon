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


# task... check if the bidirectional estimator produces better result. It would e.g. be terrific if it removed the bias in the data
E=8;
# traces = 10
# samples = 5000

runs = 10
Z = 10
T = 1
Z_min = 2
Z_max = 8

wham_list=[]
wham_bi_list = []
lr_list=[]
lr_bi_list=[]
Runs = 10
for run in range(Runs):

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
    count=0;
    while count<10:
        try:
            signal_bi, bin_centers = get_signal_bi(x_f, x_r, support, LR_bins, Z_min, Z_max, traces, samples)
            count=11
        except:
            count +=1
            pass

    ITER = 300
    slack = 1.05
    max_force = (2*pi*E*(beta**-1)/sigma)
    g_bi,z= LR_deconv_maxforce(signal_bi, 1, bin_centers,K, beta, ITER,  max_force*slack)
    g_f,z= LR_deconv_maxforce(signal_f, 1, bin_centers,K, beta, ITER, max_force*slack )
            
    g_bi = set_constant(g_bi, z, E)
    g_f = set_constant(g_f, z, E)
    lr_list.append(g_f);lr_bi_list.append(g_bi)
    
    print "deconv window: ", int(LR_bins*0.5/3.0)

    g_wham, z_valid = wham(x, wham_bins, K, len(x[0,:]), samples, support, Z_min, Z_max)
    g_wham = set_constant(g_wham, z_valid, E)
    count=0
    while count < 10:
        try:
            g_wham_bi, z_valid_bi = wham_bi(x_f, x_r, support, wham_bins, Z_min, Z_max, traces, samples, Z)
            count=11
        except:
            count += 1
            
    

    g_wham_bi = set_constant(g_wham_bi, z_valid_bi, E)

    wham_list.append(g_wham);wham_bi_list.append(g_wham_bi)
idx = arange(0, LR_bins-1, 10)
# get averages and standard deviations.
#g_wham = mean(array(wham_list), axis=0);g_wham_std =std(array(wham_list), axis=0)
#g_wham_bi =  mean(array(wham_bi_list), axis=0);g_wham_std =std(array(wham_bi_list), axis=0)  

g_f = mean(array(lr_list), axis=0);g_f_std =std(array(lr_list), axis=0)
g_bi =  mean(array(lr_bi_list), axis=0);g_bi_std =std(array(lr_bi_list), axis=0) 

ax = plt.subplot(2,1,1)
ax.set_title(r"$(U_{{0}}, s, N,N^{{LR}}_{{bins}},N^{{wham}}_{{bins}} ) = ({} [k_{{B}}T], {},{},{},{})$".format(E, samples, runs*2, LR_bins, wham_bins), size=30)
ax.errorbar(z[idx]/sigma,beta*g_f[idx],yerr=beta*g_f_std[idx]/1.0, color="blue", label=r"$g^{LR}$")
ax.errorbar(z[idx]/sigma,beta*g_bi[idx],yerr=beta*g_bi_std[idx]/1.0, color="red", label=r"$g^{LR}_{bi}$")
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
