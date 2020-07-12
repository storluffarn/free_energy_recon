#*-encoding:utf-8  -*
from numpy import *
from itertools import cycle
import matplotlib.pyplot as plt
from sys import argv
from py_utils import *
from scipy.optimize import curve_fit
from scipy.optimize import newton
from py_deconv import *

force_scale = (beta**-1)/sigma

def cutoff_function(x,x_max,cut_rate):
    """Simplest cutoff function imaginable. When x is smaller than x_max it should return 0.
    That is, in alowed regions we carry on as usual, but in forbidden regions you get no contribution"""
    for elt in x:
        if elt-x_max > 0:       # positive thing! here we want it to go to one! so that iterations slow down
            elt = 1 - e**(-cut_rate*(elt-x_max)/force_scale)
        else:
            elt = 0
    return x

def get_deriv_matrix(M,N):
    pad = ((N-M)/2) - 1
    diff = zeros((M,N))         # take our valies in N-vector, return signal shape.
    for i,row in enumerate(diff):
        row[pad + i: pad+3+i] = array([-1,0,1])
        
    return diff
    
def local_lucy(signal, extra_bins, bin_centers,K, beta, ITER,f_max):
    
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
    D = get_deriv_matrix(M,N)/(2*dz) # midpoint derivative operator.
    P_norm = get_P_sum(H,M,N)

    # cut = (N-M)/2                 # guaranteed integer. # it is also not used at all
    O = ones(N)*sum(signal)
    I = H.dot(O)
    #I_old = I                      # this is not used either
    O_old = O
    
    O = O*((H.T.dot(signal/I))/(P_norm))
    count = 0;
    while count<ITER:
        I = H.dot(O)
        O_old = O
        check = D.dot((-beta**-1)*log(O)) # here we check if the derivatives are too high.

        #check = where(check>f_max,1,0)
        check = cutoff_function(check,f_max,10)
        # we dont want to set to zero the elements that come up too big. Rather, we want to use the last signal if it was reasonable.

        
        O = O*((H.T.dot((signal/I)*(1-check)))/(P_norm))
        count+=1
        print "iter {}.".format(count)

    G = []; z = [];
    for k in range(M):
        G.append((-beta**-1)*log(O[extra_bins + k]))
        z.append(lc_bins[extra_bins+k])

    return array(G), array(z)

def get_altered_I(I, I_old, check):
    # We will use all I values that do not violate the known restrictions.
    # we shall just loop over the points in I and see what needs to be replaced with the old one
    for i,control in enumerate(check):
        if control:
            I[i] = I_old[i]
        else:
            pass
    return I
    
def local_lucy2(signal, extra_bins, bin_centers,K, beta, ITER,f_max):
    
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
    D = get_deriv_matrix(M,N)/(2*dz) # midpoint derivative operator.
    P_norm = get_P_sum(H,M,N)

    cut = (N-M)/2               # guaranteed integer.
    O = ones(N)*sum(signal)
    I = H.dot(O)
    I_old = I
    O_old = O
    
    O = O*((H.T.dot(signal/I))/(P_norm))
    count = 0;
    while count<ITER:
        O_old = O
        check = D.dot((-beta**-1)*log(O)) # here we check if the derivatives are too high.

        check = where(check>f_max,1,0) # all the elements with i should be altered
        #check = cutoff_function(check,f_max,10)
        I_old = I
        I = H.dot(O)
        I = get_altered_I(I,I_old,check)
        
        #O = O*((H.T.dot((signal/I)*(1-check)))/(P_norm))
        O = O*((H.T.dot((signal/I)))/(P_norm))

        count+=1
        print "iter {}.".format(count)

    G = []; z = [];
    for k in range(M):
        G.append((-beta**-1)*log(O[extra_bins + k]))
        z.append(lc_bins[extra_bins+k])

    return array(G), array(z)

def maxima_check(I,I_old, O, extra_bins,M, u_max):
    check = (-beta**-1)*log(O[extra_bins:extra_bins+M]) # checking in the signal range.
    check = check - check[0]                            # initializing first signal point to zero. This might not really be right though...
    check = where(check>u_max,1,0)
    for i,c in enumerate(check):
        if c:
            I[i] =I_old[i]
    return I

def local_lucy3(signal, extra_bins, bin_centers,K, beta, ITER,u_max):
    
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
    D = get_deriv_matrix(M,N)/(2*dz) # midpoint derivative operator.
    P_norm = get_P_sum(H,M,N)

    cut = (N-M)/2               # guaranteed integer.
    O = ones(N)*sum(signal)
    I = H.dot(O)
    I_old = I
    O_old = O
    
    O = O*((H.T.dot(signal/I))/(P_norm))
    count = 0;
    while count<ITER:
        O_old = O
     
        I_old = I
        I = H.dot(O)
        I = maxima_check(I,I_old, O, extra_bins, M, u_max)
        
        O = O*((H.T.dot((signal/I)))/(P_norm))
        count+=1
        print "iter {}.".format(count)

    G = []; z = [];
    for k in range(M):
        G.append((-beta**-1)*log(O[extra_bins + k]))
        z.append(lc_bins[extra_bins+k])

    return array(G), array(z)

def write_signal(signal):
    #f = open('signals/test.dat','w')
    #for item in signal:
    #      print>>signal, item
    savetxt('signals/signal.dat',signal)

def lucy_summed_max(signal, extra_bins, bin_centers,K, beta, ITER, f_m, range_tuple):
    # I think the bining is wrong, it seems to me that you get one to many bin edges, I'll just redo it...
    # Actually, new theory: need en odd number of bin centers!

    #bin_centers = delete(bin_centers,len(bin_centers)-1)
    #signal = delete(signal,len(signal)-1)
    savetxt('signals/signal.dat',signal)

    M = len(signal)
#    N = M + 2*extra_bins 
    N = M + 1

    lc_bins = zeros(N)
    dz = bin_centers[1] - bin_centers[0] # equidistant bin centers.
    bin_diff = 0.5*(bin_centers[0:len(bin_centers)-1] + bin_centers[1:len(bin_centers)])
    
    for i in range(M) :
        lc_bins[i] = bin_centers[i] - 0.5*dz
    for i in range(extra_bins) :
        lc_bins[M] = bin_centers[M-1] + 0.5*dz
    
    H = make_kernel(M,N,beta,K, bin_centers, lc_bins)
    Hsel = [0,10,19,57,150,199,500,675,750,751,752,825,998,1200,1443,1480,N-5,N-2,N-1] # for debugging
   
    # for debugging
    Hcolsum = []
    for i in range(N) :
        s = 0
        for j in range(N) :
            s += H[i,j];
        Hcolsum.append(s)
    Hrowsum = []
    for i in range(N) :
        s = 0
        for j in range(N) :
            s += H[j,i];
        Hrowsum.append(s)

    savetxt('signals/kernelcolsum.dat',Hcolsum) 
    savetxt('signals/kernelrowsum.dat',Hrowsum) 

    # for debugging
    for i in Hsel :
        Hstr = str(i)
        fname = 'signals/Hrow' + Hstr.zfill(4) + '.dat'
        savetxt(fname,H[i,:])

    P_norm = get_P_sum(H,M,N)
    #P_norm = Hrowsum

    savetxt('signals/pnorm.dat',P_norm) 

    #cut = (N-M)/2               # guaranteed integer.
    O = ones(N)*sum(signal)
    savetxt('signals/oinit.dat',O)
    
    I = H.dot(O)
    savetxt('signals/i0.dat',I)
    
    O_old = O
    n_periods = range_tuple[1] - range_tuple[0]

    signal = append(signal,0) ######## WARNING, UGLY HACK
    print shape(signal), shape(I), shape(P_norm), shape(H)
   
    O = O*((H.T.dot(signal/I))/(P_norm))
    savetxt('signals/o0.dat',-log(O))

    surface_force_max = max(diff((-beta**-1)*log(O)))/dz
    surface_force_max_old = max(diff((-beta**-1)*log(O)))/dz    # but these are just the same, why do you recalculate them?

    summed_max = summed_maxima(bin_diff, diff((-beta**-1)*log(O))/dz, range_tuple )
    summed_max_old = summed_max
    summed_f_m = (range_tuple[1] - range_tuple[0])*f_m

    max_vector = maximum_vector(bin_diff, diff((-beta**-1)*log(O))/dz,range_tuple)
    max_vector_old = max_vector
    #f_m_vector = ones(n_periods)*f_m # will not work for float inputs obviously. But in this case there is more trouble anyway.    # this is not used atm, should be commented away
    #dr = n_periods*f_m*0.05        # Add this as arg if it shows some promise. # not used either...


    # vec_size = lambda ar: sqrt(sum(ar*ar)) # just length of the vector.   # or this, damn it, comment out not used segments...
    
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
        
        count+=1
        cstr = str(count)
        fname = 'signals/oraw' + cstr.zfill(3) + '.dat'
        savetxt(fname,-log(O))
        
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
        
    
    # here. You may choose to cut even more.
    G = []; z = [];
    more_cut=0
    for k in range(N-more_cut-1):
        G.append((-beta**-1)*log(O[extra_bins + more_cut + k]))
        #G.append((-beta**-1)*log(O[extra_bins + more_cut + k]) - lc_bins[extra_bins + more_cut +k]*1E-10)       # this is some experimental slope correction!
        z.append(lc_bins[extra_bins + more_cut +k])

    return array(G), array(z)

def lucy_max(signal, extra_bins, bin_centers,K, beta, ITER, f_m, range_tuple):
    
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
#    P_norm = get_P_sum(H,M,N)

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
        
        if surface_force_max > f_m:
        #if summed_max > summed_f_m:
        #if vec_size(max_vector - f_m_vector) < dr:
            #if vec_size(max_vector - f_m_vector) < vec_size(max_vector_old-f_m_vector):
            #if abs(summed_max - summed_f_m) < abs(summed_max_old - summed_f_m):
            if abs(surface_force_max - f_m) < abs(surface_force_max_old - f_m):
                #in this case the new residual is smaller. keep just break!
                break
            else:
                O = O_old
                break
                # in this case the old one was better... make the substitution O_old ->> O and break    
        count+=1

        
        
    # here. You may choose to cut even more.
    G = []; z = [];
    more_cut=0
    for k in range(M-more_cut):
        G.append((-beta**-1)*log(O[extra_bins + more_cut + k]))
        z.append(lc_bins[extra_bins + more_cut +k])

    return array(G), array(z)


if __name__=="__main__":
    M = 10
    N = 14
    A = get_deriv_matrix(M,N)
