import numpy as np
import scipy as sp


def H0(z):
    return sp.special.hankel1(0, z, out=None)
    #return (1-1j)*np.exp(z*1j)/np.sqrt(z+0j)

def H2(z):
    return sp.special.hankel1(2, z, out=None)