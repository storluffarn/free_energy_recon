
import numpy as np

def line(scale, bins) :
    
    out = []
    
    for i in range (bins) :
        el = scale
        out.append(el)

    return out

def sinusodial(scale, periods, bins) :

    out = []

    for i in range (bins) : 
        el = np.sin(i/float(bins) * 2*np.pi * periods) + 2
        out.append(el)

    return out

## int main() I miss you <3!

bins = 1500
scale = 0.00003
periods = 6

linedata = line(scale, bins)
sinedata = sinusodial(scale,periods,bins)

np.savetxt('signals/linedata.dat',linedata)
np.savetxt('signals/sinedata.dat',sinedata)


