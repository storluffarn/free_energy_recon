#*-encoding:utf-8  -*
from numpy import *
from itertools import cycle
import matplotlib.pyplot as plt
from sys import argv
from py_utils import *
from os import system
from scipy.optimize import curve_fit

# so it scales perfectly

def f_fit(x, a,b,c):
    return a - b*abs(log(x/c))**(2.0/3.0)
def U(F_star):
    ## in units of kt
    return beta*(sigma*F_star/(2*pi))


t_list = 1.0/(array(range(1,20)))
v_list = (array(range(1,20)))
f_max = []
f = open("f_av.txt", "w");count=0;
for t in t_list:
    #system("./run -M 50 -Z 10 -E 8 -t {} -d vsV".format(t))

    system("./run -M 50 -Z 10 -E 8 -t {} -N {} -d vsV".format(t, v_list[count]))
    # second system call makes timestep k times longer if velocity is k times faster.
    
    
    x = loadtxt("vsV/x_out.dat")
    support = loadtxt("vsV/lambda_out.dat")

    f_av = average_slip_force2(x,support,K,4*10**-10)
    f_max.append(f_av)
    f.write("{} \t {} \n".format(v_list[count], f_av))
    count += 1

f.close()
popt, pcov = curve_fit(f_fit, v_list, f_max)
U_mes = U(popt[0])


# question is... what force is actually considered here....


## seems to underestimate the true value by quite a large amount. I
## the question is to which force this should apply
## lets up the number of data points but reduce the quality
## no good reason to do a check bhy hand because I do not really believe all the fits....
