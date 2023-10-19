import numpy as np
from particle_gen import *
from matrix_calc import *
from eigen import *


#Parameters
# Number of particles
N = 1000
# Exclusion radius
r_0 = 1e-4
# System's density
rho = 80.80
#--
k=1
Gamma_0 = 1
Gamma_1 = 1
#Light Direction (1 - parallel; 0 - perpendicular)
ld=0

#Generate the system of N particles on 2D space at a given density rho and an exlcusion radius of r_0 for each particle
#particle_x, particle_y = generate_normal_distribution_mk3(N, rho, r_0)
N, R, r_0 = gen_radius(N,rho,r_0)
particle_x, particle_y = generate_normal_distribution_mk2(N, R, r_0)
particles = np.column_stack((particle_x, particle_y))
plot_particles(np.sqrt(N/rho),particle_x,particle_y,N,r_0)
#Calculate the matrix with k*distance between each particle
mat = np.zeros((N,N))
for i in range(0,N):
    mat[i,i:N] = k*[np.linalg.norm(x-particles[i]) for x in particles[i:]]

#The matrix here has all values on the main diagional equals 0, and for everything else it contains the distance between each particle
mat = mat-mat.transpose()

#Construct the matrix applied to the differential equation
if ld==1:
    mat=light_par(mat,Gamma_0)
else:
    mat=light_perp_2(mat,particles,Gamma_1,N)

w,lamb=eigen_numpy(mat)
print(np.max(w))
cen_mass=sum(particles)/N
k=[np.linalg.norm(x-cen_mass) for x in particles]

plt.scatter(k,np.abs(w)**2)
plt.show()





