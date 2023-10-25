import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
import time
import random

#The whole system will match the desirable density, but it could have regions where the density is different
#Here the maximum possible density will be 1/(pi*r²_0)
def gen_radius(N,rho,r_0):
    #rho - Particles density of the system
    area = N/rho
    if N*(np.pi*(r_0**2)) >= area:
        raise ValueError("Exclusion area is bigger than avaliabe space to include all particles.")
    R = np.sqrt(area/np.pi)
    return N,R,r_0

#Version mk1
#Given a radius R, the system choices random values for x and y within the maximum radius
#After that the tree checks if the position is already within the exclusion range of another particle AND if it's inside the area requested
def generate_normal_distribution_mk1(N, R, r_0):
    #N - Number of particles on the system
    #r_0 - Individual exclusion radius
    #R - Radius of area where the partciles should be placed   
    
    #Creates coordinates x and y
    particle_x = []
    particle_y = []
    
    #Tree contains all the particles inside, which are gonna be used for check if those new candidates are within exlusion radius or not
    tree = None

    #Check if the number of particles match the desirable N
    while len(particle_x) < N:
        x = np.random.uniform(-R, R)
        y = np.random.uniform(-R, R)

        #Check first if the candidate coordinates are whitin exclusion radius for any particle inside the tree, then check if it's inside the requested area
        if (tree is None or not any(tree.query_ball_point([x, y], r_0))) and (x**2 + y**2 <= R**2):
            particle_x.append(x)
            particle_y.append(y)
            tree = cKDTree(list(zip(particle_x, particle_y)))
    return particle_x, particle_y

#Version mk2
#Now the system choices based on polar coordinates, where the analysis to check if the position is inside the range can be eliminated
def generate_normal_distribution_mk2(N, R, r_0):
    #N - Number of particles on the system
    #r_0 - Individual exclusion radius
    #R - Radius of area where the partciles should be placed   
    #Creates coordinates x and y
    particle_x = []
    particle_y = []
    #Tree contains all the particles inside, which are gonna be used for check if those new candidates are within exlusion radius or not
    tree = None
    #Check if the number of particles match the desirable N
    while len(particle_x) < N:
        r = np.sqrt(np.random.uniform(0, R**2))
        theta = np.random.uniform(0, 2*np.pi)
        x,y = r*np.cos(theta),r*np.sin(theta)
        #Check first if the candidate coordinates are whitin exclusion radius for any particle inside the tree
        if (tree is None or not any(tree.query_ball_point([x, y], r_0))):
            particle_x.append(x)
            particle_y.append(y)
            tree = cKDTree(list(zip(particle_x, particle_y)))
    return particle_x, particle_y

#This way instead of a circular area, there's a square where each side is L, in such way that N/L² = rho
#Now the area is divided in small grids, each one with sides equals the exclusion radius
#The system choices for x and y some index related to one of these grids. After the choice, this option is deleted,
#so there shouldn't be any other candidate trying to get the grid already occupied

def generate_normal_distribution_mk3(N, rho, r_0):
    L = np.sqrt(N/rho)
    #Grants at least a minium value to create the grids
    if r_0==0:
        r_0=1e-3
    n_grids= int(np.sqrt(L/r_0))
    if n_grids**2 < N:
        raise ValueError("Exclusion radius r_0 must be smaller or density should be reduced.")
    #Create the index for all grids
    pool=[*range(n_grids**2)]
    particle_x,particle_y = [],[]
    
    while len(particle_x) < N:
        x = random.choice(pool)
        pool.remove(x)    
        particle_x.append(L*(x%n_grids)/n_grids)
        particle_y.append(L*(x//n_grids)/n_grids)
        
    return particle_x,particle_y

# Plot the particles
def plot_particles(R,particle_x,particle_y,N,r_0):
    plt.figure(figsize=(8, 8))
    plt.scatter(particle_x,particle_y, s=10)
    plt.xlim(-R, R)
    plt.ylim(-R, R)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title(f'N={N}, R={R}, r_0={r_0}')
    plt.xlabel('X-coordinate')
    plt.ylabel('Y-coordinate')
    plt.grid()
    plt.show()

