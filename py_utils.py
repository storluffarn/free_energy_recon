#*-encoding:utf-8  -*
from numpy import *
from itertools import cycle
import matplotlib.pyplot as plt
from sys import argv
kb = 1.38064852*10**-23
beta = (300*kb)**-1
sigma = 0.3*10**(-9)
K = 2

## First order of business is to find the average slip force for a given trajectory.

def averge_slip_force(x,support,K,thresh):
    """Given a coordinate trajectory, returns the average value of the slipforce."""
    force_trace = x*0
    for k in range(len(x[0,:])):
        force_trace[:,k] = -K*( x[:,k]-support)

    ## We should be cognizant of the fact that this need to be tested. Make a simplistic "if you move by this amount" trial and save the points where it triggers.
    s_list = []
    x_list = []

    # This is important, this is the slipping criterion! 
    for k in range(len(x[0,:])):
        for j in range(1,len(x[:,0])):
            if -force_trace[j,k] + force_trace[j-1,k] > thresh:
                ## In this case the force was significantly reduced, so j-1 value ought to be where it slipped.
                x_list.append(force_trace[j-1,k])
                s_list.append(support[j-1])

    return array(s_list), array(x_list)

def average_slip_force2(x,support,K,thresh):
    s_m,f_m = averge_slip_force(x,support,K,thresh) # what is the syntactic meaning of comma here?
    return mean(f_m)



def cluster_to_min(cluster):
    """cluster is a list of lists, clusters. """
    min_list = []
    for l in cluster:
        min_list.append(min(array(l)))
    return array(min_list)

def first_in_cluster(s, x, MAX_SPACING ):
    """Given an array of forces before slipping and an array of corresponding support-positions, returns a list of the first such occurence in every slip event"""

    cluster_list = []           # list of lists.
    cluster_list.append([])

    cluster_counter = 0;
    last_pos = s[0]
    cluster_list[cluster_counter].append(last_pos)

    for i in range(1, len(x)):
        if (s[i]-last_pos) > MAX_SPACING:
            cluster_list.append([]); cluster_counter+=1;
            cluster_list[cluster_counter].append(x[i])
            last_pos = s[i]
            #cluster_list[cluster_counter].append(last_pos)
            
        else:
            cluster_list[cluster_counter].append(x[i])  # why append rather than 
            # updating?
            last_pos = s[i]
            #cluster_list[cluster_counter].append(last_pos)
            
    return cluster_list

def get_first_slip_forces(s,x, MAX_SPACING):

    s_list = []
    x_list = []
    last_pos = s[0]
    s_list.append(s[0]);x_list.append(x[0])

    for i in range(1,len(x)):
        if (s[i]-s[i-1]) > MAX_SPACING:
            s_list.append(s[i]); x_list.append(x[i])
        else:
            pass
            # if not, do nothing and just move on
    return array(s_list), array(x_list)

if __name__ == "__main__":  # this sorts of run the program then? how does __main__ work in python?
    x = loadtxt("vsV/x_out.dat")
    support = loadtxt("vsV/lambda_out.dat")
    s,f_m = averge_slip_force(x,support,K,4*10**-10)
    # cluster = first_in_cluster(s, f_m, sigma*1.5 )
    # min_array = cluster_to_min(cluster) # this is wrong though... you do not need all of this shit.

    #s_m, f_m = get_first_slip_forces(s,f_m, sigma*2.5)

    F = averge_slip_force2(x,support,K,4*10**-10)
    # ax = plt.subplot(1,1,1)
    # for k in range(len(x[0,:])):
    #     ax.plot(support, -K*(x[:,k]-support), color="red", alpha=0.4)
    # ax.plot(s,f_m, "o", ms=14)
