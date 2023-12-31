import numpy as np
from particle_gen import *
from matrix_calc import *
from eigen import *


#Parameters
# Number of particles
N = 10000
# Exclusion radius
r_0 = 0
# System's density
rho = 1
#--
k=1
Gamma_0 = 1
Gamma_1 = 1
#Light Direction (1 - parallel; 0 - perpendicular)
ld=1

#Generate the system of N particles on 2D space at a given density rho and an exlcusion radius of r_0 for each particle
#particle_x, particle_y = generate_normal_distribution_mk2(N, R, r_0)
N, R, r_0 = gen_radius(N,rho,r_0)
particle_x, particle_y = generate_normal_distribution_mk2(N, R, r_0)

particles = np.column_stack((particle_x, particle_y))
#plot_particles(np.sqrt(N/rho),particle_x,particle_y,N,r_0)

#Calculate the matrix with k*distance between each particle
mat = np.zeros((N,N))
for i in range(0,N):
    mat[i,i:N] = k*[np.linalg.norm(x-particles[i]) for x in particles[i:]]
#The matrix here has all values on the main diagional equals 0, and for everything else it contains the distance between each particle
mat = mat+mat.transpose()

#Construct the matrix applied to the differential equation

mat=light_par(mat,Gamma_0)

lamb,w=eigen_numpy(mat)
print(lamb)
plt.scatter(lamb.imag,lamb.real)
plt.yscale('log')
plt.show()




