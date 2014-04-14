from math import log
import numpy as np
import matplotlib.pyplot as plt


def T_binomial(p, g, m):
    alpha = 5.000000000000001e-07
    beta = 9.313225746154785e-12
    return log(p,2)*(alpha + beta)



def T_scatter_recd_allgather(p, g, m):
    alpha = 5.000000000000001e-07
    beta = 9.313225746154785e-12
    if (g==1 or g==p):
        return (2*log(p,2))*alpha + 2*m*(1-1.0/p)*beta
    else:
        return ( 2*log(p,2))*alpha + 2*m*( 2 - 1.0/g - g/p )*beta


def T_scatter_ring_allgather(p, g, m):
    alpha = 5.000000000000001e-07
    beta = 9.313225746154785e-12
    if (g==1 or g==p):
        return (log(p,2) + p -1)*alpha + 2*m*(1-1.0/p)*beta
    else:
        return ( log(p,2) + g + p/g -2 )*alpha + 2*m*( 2 - 1.0/g - g/p )*beta



powers = range(0,15)

groups=np.power(len(powers)*[2], powers)


plt.plot(groups, [T_scatter_ring_allgather(16384, g, 16*1024) for g in groups])

plt.xscale('log', basex=2)

plt.show()

