import numpy as np

def H0(z):
    return (1-1j)*np.exp(z*1j)/np.sqrt(z+0j)

def H2(z):
    return -(1-1j)*np.exp(z*1j)/np.sqrt(z+0j)