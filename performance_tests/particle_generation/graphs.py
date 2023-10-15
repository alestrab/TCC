import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit

cwd = os.getcwd()
print(cwd)
files = [f for f in os.listdir(cwd+'\data')]
i=0

def fit(x,a,b,c):
    return a*x**3+b*x**2+c*x

for f in files:
    if i==0:
        color='red'
        label = 'mk1 - Cartesian Coordinates + Tree Confirmation'
    elif i==1:
        color='blue'
        label = 'mk2 - Polar Coordinates + Tree Confirmation'
    else:
        label = 'mk3 - Small Grids'
        color='green'
    i+=1
    data = np.loadtxt(cwd+'\data\\'+f)
    x,y = data[:,0],data[:,1]
    plt.scatter(x,y,color=color,label=label)
    xdata = np.linspace(0, 50000, 10000)
    popt, pcov = curve_fit(fit, x, y)
    plt.plot(xdata,fit(xdata,*popt),color=color)
plt.xlabel('Number of Particles',fontsize=15)
plt.ylabel('Time to Generate the System [s]',fontsize=15)
plt.legend()
plt.grid()
plt.show()

