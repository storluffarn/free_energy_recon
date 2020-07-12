#*-encoding:utf-8  -*
from numpy import *
from itertools import cycle
import matplotlib.pyplot as plt
from sys import argv
from py_utils import *
from scipy.optimize import curve_fit
from scipy.optimize import newton

# som first and foremoste, you need to create a signal for the algorithm to take.
# The problem is in the binning process. T
# There is something wrong with the WHAM case. It does not give a good fit for the low energies, where it is supposed to be excellent.


def get_signal(x, support, bins, z_min, z_max):

    ## There is an error here that causes regular, fast oscillations. Find it to remove this taunting message!
    
    z_min = sigma*z_min
    z_max = sigma*z_max
    d_bin = (z_max-z_min)/bins
    
    traces = len(x[0,:])
    samples = len(x[:,0])

    print "traces,samples", traces,samples
    bin_edges = zeros(bins+1); bin_centers = zeros(bins)
    for i in range(bins+1):
        bin_edges[i] = z_min + d_bin*i
    for i in range(bins):
        bin_centers[i] = 0.5*(bin_edges[i] + bin_edges[i+1])
    
    ext_work = zeros((traces, samples))
    exp_average = zeros(samples)
    signal = zeros(bins)
    signal_counter = signal*0

    
    for k in range(traces):
        ext_work[k,0] = 0;      # set first to zero
        for i in range(1,samples):
            ext_work[k,i] = ext_work[k,i-1] - 0.5*K*( (x[i,k]-support[i]) + (x[i-1,k]- support[i-1]) )*(support[i] - support[i-1])
    for i in range(samples):
        for k in range(traces):
            exp_average[i] += exp(-beta*ext_work[k,i])
        exp_average[i] = exp_average[i]/traces
    for i in range(samples):
        for l in range(bins):
            if bin_edges[l] < support[i] and bin_edges[l+1] >= support[i]:
                signal[l] += exp_average[i]
                signal_counter[l] += 1.0
    for l in range(bins):
        if signal_counter[l] != 0:
            signal[l] = signal[l]/signal_counter[l]

    ## not normalized now...
    return signal, bin_centers

def get_signal_for_all_samples(x, support, z_min, z_max):
    """Instead of histogramming arbitrary """
    traces = len(x[0,:])
    samples = len(x[:,0])
    x_new = zeros((samples,traces))
    
    bins = array([s for s in support if s > z_min*sigma and s <= z_max*sigma])
    for k in range(traces):
        x_new[:,k] = array([x[i,k] for i in range(samples) if support[i] > z_min*sigma and support[i] <= z_max*sigma])

    return bins, x_new
    
def surface_potential(bin_centers, E):
    return (E*beta**-1)*(1-cos(bin_centers*2*pi/sigma))
# def surface_potential(bin_centers, E):

#     return E*(64*(beta**-1)*(sigma**-4))*(  0.25*bin_centers**(4) - 0.5*sigma*bin_centers**(3) + 0.25*(sigma**2)*(bin_centers**2)   )

def summed_maxima(x,z,range_tuple):
    """Given a reconstructed surface, returns the sum of the maxima in each period. Range tuple is of form (z_min, z_max) in units of sigma"""
    n_periods = range_tuple[1] - range_tuple[0]
    S = 0;
    for i in range(n_periods):
        t = where(logical_and(greater_equal(x, sigma*(i + range_tuple[0])), less_equal(x, sigma*(i + 1 + range_tuple[0] )) ))
        S +=max(z[t])
    return S
def maximum_vector(x,z,range_tuple):
    n_periods = range_tuple[1] - range_tuple[0]
    S = [];
    for i in range(n_periods):
        t = where(logical_and(greater_equal(x, sigma*(i + range_tuple[0])), less_equal(x, sigma*(i + 1 + range_tuple[0] )) ))
        S.append(max(z[t]))
    return array(S)
    
def make_pbc_kernel(beta, bin_centers,N, extra_bins, d_bin):

    # if 
    kernel = zeros((N,N))
    norm_vector = zeros(N)
    for i in range(N):
        # create the first element first
        kernel[i,i] = 1;        # because distance is zero
        for j in range(1,extra_bins):
            kernel[i, (i-j)%N] += exp(-0.5*beta*K*((j*d_bin)**2))
            kernel[i, (i+j)%N] += exp(-0.5*beta*K*((j*d_bin)**2))
    return kernel
            
def make_kernel(M,N,beta,K, bin_centers, lc_bins):
    kernel = zeros((M,N))
    norm_vector = zeros(N)
    for i in range(M):
        for j in range(N):
            kernel[i,j] = exp(-0.5*beta*K*( (bin_centers[i] - lc_bins[j])**2 ))
            norm_vector[j] += exp(-0.5*beta*K*( (bin_centers[i] - lc_bins[j])**2 ))
    for i in range(M):
        for j in range(N):
            kernel[i,j] = kernel[i,j]/norm_vector[j]
    return kernel
def get_P_sum(H,M,N):
    ## this really should make a different for a normalized kernel.... but.. im not confident enough
    ## to make this call atm.
    P = zeros(N)
    for j in range(N):
        for i in range(M):
            P[j] += H[i,j]
    return P

def set_constant(g, z,E):
    const = -g[0] + surface_potential(z[0],E)
    return g + const
    
def LR_deconv(signal, extra_bins, bin_centers,K, beta, ITER):
    #ITER_MAX=20
    M = len(signal)
    N = M + 2*extra_bins
    lc_bins = zeros(N)
    dz = bin_centers[1] - bin_centers[0] # equidistant bin centers.

    for i in range(M):
        lc_bins[extra_bins + i] = bin_centers[i]
    for i in range(extra_bins):
        lc_bins[M + i + extra_bins] = bin_centers[M-1] + (i+1)*dz
        lc_bins[extra_bins - i - 1] = bin_centers[0] - (i+1)*dz
    H = make_kernel(M,N,beta,K, bin_centers, lc_bins)
    P_norm = get_P_sum(H,M,N)

    cut = (N-M)/2               # guaranteed integer.
    O = ones(N)*sum(signal)
    I = H.dot(O)
    O_old = O

    
    O = O*((H.T.dot(signal/I))/(P_norm))
    count = 0;
    while count<ITER:
        
        I = H.dot(O)            # Data guess, residuals I-signal shoud be small.
        O_old = O
        #O = O*(H.T.dot(signal/I))
        O = O*((H.T.dot(signal/I))/(P_norm))
        count+=1


    G = []; z = [];
    for k in range(M):
        G.append((-beta**-1)*log(O[extra_bins + k]))
        z.append(lc_bins[extra_bins+k])

    return array(G), array(z)
def LR_deconv_maxforce(signal, extra_bins, bin_centers,K, beta, ITER, f_m, range_tuple):
    #ITER_MAX=20
    # 3 different strategies here makes this code pretty slow.
    # problem with the residual method is that its very possible to completely miss the minima.
    M = len(signal)
    N = M + 2*extra_bins
    lc_bins = zeros(N)
    dz = bin_centers[1] - bin_centers[0] # equidistant bin centers.
    bin_diff = 0.5*(bin_centers[0:len(bin_centers)-1] + bin_centers[1:len(bin_centers)])
    
    for i in range(M):
        lc_bins[extra_bins + i] = bin_centers[i]
    for i in range(extra_bins):
        lc_bins[M + i + extra_bins] = bin_centers[M-1] + (i+1)*dz
        lc_bins[extra_bins - i - 1] = bin_centers[0] - (i+1)*dz
    H = make_kernel(M,N,beta,K, bin_centers, lc_bins)
    P_norm = get_P_sum(H,M,N)

    cut = (N-M)/2               # guaranteed integer.
    O = ones(N)*sum(signal)
    I = H.dot(O)
    O_old = O
    n_periods = range_tuple[1] - range_tuple[0]
    
    O = O*((H.T.dot(signal/I))/(P_norm))

    surface_force_max = max(diff((-beta**-1)*log(O)))/dz
    surface_force_max_old = max(diff((-beta**-1)*log(O)))/dz

    summed_max = summed_maxima(bin_diff, diff((-beta**-1)*log(O))/dz, range_tuple )
    summed_max_old = summed_max
    summed_f_m = (range_tuple[1] - range_tuple[0])*f_m

    max_vector = maximum_vector(bin_diff, diff((-beta**-1)*log(O))/dz,range_tuple)
    max_vector_old = max_vector
    f_m_vector = ones(n_periods)*f_m # will not work for float inputs obviously. But in this case there is more trouble anyway.
    dr = n_periods*f_m*0.05                # Add this as arg if it shows some promise.


    vec_size = lambda ar: sqrt(sum(ar*ar)) # just length of the vector.
    count = 0;
    while count<ITER:
        
        I = H.dot(O)            # Data guess, residuals I-signal shoud be small.
        O_old = O
        #O = O*(H.T.dot(signal/I))
        O = O*((H.T.dot(signal/I))/(P_norm))

        # Max criterion
        surface_force_max_old = surface_force_max
        surface_force_max = max(diff((-beta**-1)*log(O)))/dz

        # Summed Max Criterion.
        summed_max_old = summed_max
        summed_max = summed_maxima(bin_diff, diff((-beta**-1)*log(O))/dz, range_tuple )

        # Residual Criterion.
        max_vector_old = max_vector
        max_vector = maximum_vector(bin_diff, diff((-beta**-1)*log(O))/dz,range_tuple)
        
        #if surface_force_max > f_m:
        if summed_max > summed_f_m:
        #if vec_size(max_vector - f_m_vector) < dr:
            #if vec_size(max_vector - f_m_vector) < vec_size(max_vector_old-f_m_vector):
            if abs(summed_max - summed_f_m) < abs(summed_max_old - summed_f_m):
            #if abs(surface_force_max - f_m) < abs(surface_force_max_old - f_m):
                #in this case the new residual is smaller. keep just break!
                break
            else:
                O = O_old
                break
                # in this case the old one was better... make the substitution O_old ->> O and break
        
        print "max surface force of reconstructed potential = ", surface_force_max, "({})".format(f_m)
        print "Summed max = {}, target = {}".format(summed_max, summed_f_m)
        print "Residual: {}, target: {}".format(vec_size(max_vector - f_m_vector), dr)
        count+=1

        
        
    # here. You may choose to cut even more.
    G = []; z = [];
    more_cut=0
    for k in range(M-more_cut):
        G.append((-beta**-1)*log(O[extra_bins + more_cut + k]))
        z.append(lc_bins[extra_bins + more_cut +k])

    return array(G), array(z)

def LR_deconv_maxforce_pbc(signal, extra_bins, bin_centers,K, beta, ITER, f_m):
    #ITER_MAX=20
    N = len(signal)
    dz = bin_centers[1] - bin_centers[0] # equidistant bin centers.

    H = make_pbc_kernel(beta, bin_centers, N, extra_bins, dz)
    #H = make_kernel(M,N,beta,K, bin_centers, lc_bins)
    P_norm = get_P_sum(H,N,N)

    O = ones(N)*sum(signal)
    I = H.dot(O)
    O_old = O

    
    O = O*((H.T.dot(signal/I))/(P_norm))

    surface_force_max = max(diff((-beta**-1)*log(O)))/dz
    surface_force_max_old = max(diff((-beta**-1)*log(O)))/dz
    count = 0;
    while count<ITER:
        
        I = H.dot(O)            # Data guess, residuals I-signal shoud be small.
        O_old = O
        #O = O*(H.T.dot(signal/I))
        O = O*((H.T.dot(signal/I))/(P_norm))
        surface_force_max_old = surface_force_max
        surface_force_max = max(diff((-beta**-1)*log(O)))/dz

        # uncomment if you want global maximum force stop. This must be made local in some way, in a systematic fashion.
        if surface_force_max > f_m:
            # if the maximum surface force has become stronger, we should check which of the iterations comes closest to the real thing and return that
            if abs(surface_force_max - f_m) < abs(surface_force_max_old - f_m):
                #in this case the new residual is smaller. keep just break!
                break
            else:
                O = O_old
                break
                # in this case the old one was better... make the substitution O_old ->> O and break
        
        # print "max surface force of reconstructed potential = ", surface_force_max, "({})".format(f_m)
        count+=1

        
    G = []; z = [];
    for k in range(N):
        G.append((-beta**-1)*log(O[k]))
        z.append(bin_centers[k])

    return array(G), array(z)



def self_consistant_eq(dF, w_f, w_r, traces):
    top = 0
    bottom = 0;
    for k in range(traces):
        top += (1 + exp((w_f[k] - dF)))**(-1)
        bottom += (1 + exp((w_r[k] + dF)))**(-1)
    return 1 - (top/bottom)
        
def BAR_dF(ext_work_f, ext_work_r, traces, samples):
    """Implements the crooks version of Bennet Accept Ratio (BAR) to find a more realistic value for the dF,
    (hopefully this will be as close as possible to dF)"""
    w_f = ext_work_f[:, samples-1]*beta # correct units! otherwise it will NOT converge at all
    w_r = ext_work_r[:, samples-1]*beta

    temp = newton(self_consistant_eq,0, args=(w_f, w_r, traces))
    print "Residual from Bar fit: ", self_consistant_eq(temp, w_f, w_r, traces) 
    return (beta**-1)*temp    
    # return newton(0,self_consistant_eq, args=(w_f, w_r, traces))
    
    
def jarz_dF(ext_work_f, ext_work_r, traces,samples):
    """If the process is entirely symmetric, it seems unnecessary to calculate the free energy change self-consistantly"""
    temp_exp = 0
    for k in range(traces):
        temp_exp += exp(-beta*ext_work_f[k, samples-1])
        temp_exp += exp(-beta*ext_work_r[k, samples-1])
    temp_exp = temp_exp/(2.0*traces)

    return -(beta**-1)*log(temp_exp)
    
    
def get_signal_bi(x_f, x_r, support, bins, z_min, z_max, traces, samples):
    """ In this case there will be the same amoujnt of trajectories in either direction."""
    # the work W(tau-t --> tau)_r can be calculated with a simple difference. We probably do not need alof of new structure to make this work
    # Alot of inefficiencies are in the code, mostly because it makes it easier to change if I want a amore flexible code.
    
    z_min = sigma*z_min
    z_max = sigma*z_max
    d_bin = (z_max-z_min)/bins

    bin_edges = zeros(bins+1); bin_centers = zeros(bins)
    for i in range(bins+1):
        bin_edges[i] = z_min + d_bin*i
    for i in range(bins):
        bin_centers[i] = 0.5*(bin_edges[i] + bin_edges[i+1])

    # strategy: calculate the signal at every smaple point and THEN average them into bins...
    ext_work_f = zeros((traces, samples))
    exp_work_f = zeros(samples)
    ext_work_r = ext_work_f*0; exp_work_r = exp_work_f*0
    exp_work_tot = exp_work_f*0

    signal = zeros(bins)
    signal_counter = zeros(bins)

        
    for k in range(traces):
        ext_work_f[k,0] = 0;      # set first to zero
        ext_work_r[k,0] = 0;
        for i in range(1,samples):
             ext_work_f[k,i] = ext_work_f[k,i-1] - 0.5*K*( (x_f[i,k]-support[i]) + (x_f[i-1,k]- support[i-1]) )*(support[i] - support[i-1])
             ext_work_r[k,i] = ext_work_r[k,i-1] - 0.5*K*( (x_r[i,k]-support[i]) + (x_r[i-1,k]- support[i-1]) )*(support[i] - support[i-1])
        # because the works are really the same here, because problem is symmetric, the combination will be different though.

    #dF = 0;                     # this can be relaxed later.
    # dF = jarz_dF(ext_work_f, ext_work_r, traces,samples)
    # print "jarz dF: ", dF
    dF = BAR_dF(ext_work_f, ext_work_r, traces, samples)
    print "BAR dF: ", dF
    
    for i in range(samples):
        for k in range(traces):
            # the exponential averages appearing in Eq 7 of Minhs article. So, the standard for forward dynamics and
            # exp(beta*W(tau-t ---> tau)_r) for 
            exp_work_f[i] += traces*exp(-beta*ext_work_f[k,i])/(traces + traces*exp(-beta*(ext_work_f[k, samples-1] - dF))) 
            temp_work = ext_work_r[k, samples-1] - ext_work_r[k, samples - 1 - i] # so that when i = samples last one says ext_work_r[k, samples - 1 - (samples-1)] = ext_work_r[0]
            ## the work W(tau-t ----> tau) = W(0 --> tau) - W(0 --> tau-t), should be
            exp_work_r[i] += traces*exp(beta*temp_work)/(traces + traces*(exp(beta*(ext_work_r[k, samples-1] + dF))))

        exp_work_tot[i] = (exp_work_f[i]/traces) + (exp_work_r[i]/traces)

    for i in range(samples):
        for l in range(bins):
            if bin_edges[l] < support[i] and bin_edges[l+1] > support[i]:
                signal[l] += exp_work_tot[i]
                signal_counter[l] += 1
    for l in range(bins):
        if signal_counter[l] != 0:
            signal[l] = signal[l]/signal_counter[l]

    
    # return exp_work_tot
    return signal, bin_centers

def remove_zeros(G_exp, bin_centers):
    G = []
    z = []
    for l in range(len(G_exp)):
        if G_exp[l] != 0:
            G.append(-(beta**-1)*log(G_exp[l]))
            z.append(bin_centers[l])
    return array(G), array(z)
def wham(x, bins, K, traces, samples, support, z_min, z_max):
    """ Implements a classic weighted histogram analysis method based on forward trajectories """
    # Shit is incorrectly normalized!

    z_min = sigma*z_min
    z_max = sigma*z_max
    d_bin = (z_max-z_min)/bins

    bin_edges = zeros(bins+1); bin_centers = zeros(bins)
    for i in range(bins+1):
        bin_edges[i] = z_min + d_bin*i
    for i in range(bins):
        bin_centers[i] = 0.5*(bin_edges[i] + bin_edges[i+1])
    
    hist = zeros((bins, samples))
    hist_counter = hist*0
    weights = zeros(samples)
    nominator = zeros(bins)
    denominator = nominator*0

    ext_work = zeros((traces, samples))
    
    for k in range(traces):
        ext_work[k,0] = 0;      # set first to zero
        for i in range(1,samples):
            ext_work[k,i] = ext_work[k,i-1] - 0.5*K*( (x[i,k]-support[i]) + (x[i-1,k]- support[i-1]) )*(support[i] - support[i-1])

    
    for i in range(samples):
        for k in range(traces):
            weights[i] += exp(-beta*ext_work[k,i])
            for l in range(bins):
                if (x[i,k] > bin_edges[l] and x[i,k] <= bin_edges[l+1]):
                    #print "Hit!"
                    hist[l,i] += exp(-beta*ext_work[k,i])/traces # normalizing on every term, for some reason
                    #hist_counter[l,i] += 1.0
        # loop over traces done, normalizing where there are hits, and normalizing all weights (guaranteed hits).
        #for l in range(bins):
            #if hist_counter[l,i] != 0:
                #hist[l,i] = hist[l,i]/hist_counter[l,i]
        weights[i] = weights[i]/traces
    
    V = lambda z,l: exp(-beta*0.5*K*((z-l)**2))
    G_exp = zeros(bins)
    for l in range(bins):
        for i in range(samples):
            nominator[l] += hist[l,i]/weights[i] 
            denominator[l] += V(bin_centers[l], support[i])/weights[i]
        G_exp[l] = nominator[l]/denominator[l]
    
    G, z = remove_zeros(G_exp, bin_centers)
    #print "Reconstructing ({},{}) using {} bins ".format(min(x.flatten())*10**9, max(x.flatten())*10**9, bins)
    return G,z
def wham_bi(x_f, x_r, support, bins, z_min, z_max, traces, samples, Z):
    """ In this case there will be the same amoujnt of trajectories in either direction."""
    # the work W(tau-t --> tau)_r can be calculated with a simple difference. We probably do not need alof of new structure to make this work
    # Alot of inefficiencies are in the code, mostly because it makes it easier to change if I want a amore flexible code.
    
    z_min = sigma*z_min
    z_max = sigma*z_max
    d_bin = (z_max-z_min)/bins

    bin_edges = zeros(bins+1); bin_centers = zeros(bins)
    for i in range(bins+1):
        bin_edges[i] = z_min + d_bin*i
    for i in range(bins):
        bin_centers[i] = 0.5*(bin_edges[i] + bin_edges[i+1])

    # strategy: calculate the signal at every smaple point and THEN average them into bins...
    ext_work_f = zeros((traces, samples))
    exp_work_f = zeros(samples)
    ext_work_r = ext_work_f*0; exp_work_r = exp_work_f*0
    exp_work_tot = exp_work_f*0

    signal = zeros(bins)
    signal_counter = zeros(bins)

        
    for k in range(traces):
        ext_work_f[k,0] = 0;      # set first to zero
        ext_work_r[k,0] = 0;
        for i in range(1,samples):
             ext_work_f[k,i] = ext_work_f[k,i-1] - 0.5*K*( (x_f[i,k]-support[i]) + (x_f[i-1,k]- support[i-1]) )*(support[i] - support[i-1])
             ext_work_r[k,i] = ext_work_r[k,i-1] - 0.5*K*( (x_r[i,k]-support[i]) + (x_r[i-1,k]- support[i-1]) )*(support[i] - support[i-1])
        # because the works are really the same here, because problem is symmetric, the combination will be different though.

    #dF = 0;                     # this can be relaxed later.
    # dF = jarz_dF(ext_work_f, ext_work_r, traces,samples)
    # print "jarz dF: ", dF
    dF = BAR_dF(ext_work_f, ext_work_r, traces, samples)
    print "BAR dF: ", dF
    
    for i in range(samples):
        for k in range(traces):
            # the exponential averages appearing in Eq 7 of Minhs article. So, the standard for forward dynamics and
            # exp(beta*W(tau-t ---> tau)_r) for 
            exp_work_f[i] += traces*exp(-beta*ext_work_f[k,i])/(traces + traces*exp(-beta*(ext_work_f[k, samples-1] - dF))) 
            temp_work = ext_work_r[k, samples-1] - ext_work_r[k, samples - 1 - i] # so that when i = samples last one says ext_work_r[k, samples - 1 - (samples-1)] = ext_work_r[0]
            ## the work W(tau-t ----> tau) = W(0 --> tau) - W(0 --> tau-t), should be
            exp_work_r[i] += traces*exp(beta*temp_work)/(traces + traces*(exp(beta*(ext_work_r[k, samples-1] + dF))))

        exp_work_tot[i] = (exp_work_f[i]/traces) + (exp_work_r[i]/traces)

    hist_f = zeros((bins, samples)); hist_r = hist_f*0
    hist_counter_f = hist_f*0; hist_counter_r = hist_f*0

    #nominator_f = zeros(bins); nominator_r = nominator_f*0;
    nominator = zeros(bins)
    denom = nominator*0
    g_exp = nominator*0
    # All this is probbly not necessary, but it is a more general way and thus easier to change to REAL reverse dynamics.


    # shouldnt you jsuts reverse it instead?
    x_r_mirror = Z*sigma - x_r  # Z is maximum support position in units of sigma.
    
    for i in range(samples):
        for k in range(traces):
            for l in range(bins):
                if (x_f[i,k] > bin_edges[l] and x_f[i,k] <= bin_edges[l+1]):
                    hist_f[l,i] += (traces*exp(-beta*ext_work_f[k,i])/(traces + traces*exp(-beta*(ext_work_f[k, samples-1]-dF))))/traces # normalizing on every term, for some reason
                    #hist_counter_f[l,i] += 1.0

                # this is trickier than it loooks...
                if (x_r_mirror[samples -1 -i,k] > bin_edges[l] and x_r_mirror[samples -1 -i ,k] <= bin_edges[l+1]):
                    temp_work = ext_work_r[k, samples-1] - ext_work_r[k, samples - 1 - i]
                    hist_r[l,i] +=  (traces*exp(beta*temp_work)/(traces + traces*exp(beta*(ext_work_r[k,samples-1]+dF))))/traces# normalizing on every term, for some reason
                    #hist_counter_r[l,i] += 1.0
                
                    
        # for l in range(bins):
        #     if hist_counter_f[l,i] != 0:
        #         hist_f[l,i] = hist_f[l,i]/hist_counter_f[l,i]
        #     if hist_counter_r[l,i] != 0:    
        #         hist_r[l,i] = hist_r[l,i]/hist_counter_r[l,i]
        
        # now we have to just sum over all time and put it together
    V = lambda z,l: exp(-beta*0.5*K*((z-l)**2))
    for l in range(bins):
        for i in range(samples):
            nominator[l] += (hist_f[l,i] + hist_r[l,i])*(exp_work_tot[i]**-1)
            denom[l] += V(bin_centers[l], support[i])*(exp_work_tot[i]**-1)
        g_exp[l] = nominator[l]/denom[l]

    g,z = remove_zeros(g_exp, bin_centers)
        
    return g,z
