from numpy import fill_diagonal,zeros,ones
from hankel_func import *

def light_par(mat,Gamma_0):
    #Perform Hankel function for all the matrix
    #Add 1 to the main diagonal just to don't crash the vectorized calculus
    fill_diagonal(mat, 1)
    mat = H0(mat)
    #Correct the value for the main diagonal
    fill_diagonal(mat, 1)
    #mat = mat*(-Gamma_0/2)
    return mat



def light_perp(mat,particles,Gamma_1,N):
    #Perform Hankel function for all the matrix
    #Add 1 to the main diagonal just to don't crash the vectorized calculus
    fill_diagonal(mat, 1)
    
    mat_hank_zero=H0(mat)
    mat_hank_two = H2(mat)
    
    #Create matrix with coordinates distance and sum the components
    mat_dist = zeros((N,N))
    #Here phi_jl = (y[j]-y[l])/(x[j]-x[l]), and so, phi_jl=phi_lj
    #Also, the line 'i' must contain phi_i0,phi_i1, ..., phi_iN
    #So, on line 'i' the particle 'i' coordinates remain fixed, while the others vary for all particles in the system
    for i in range(0,N):
        mat_dist[i,i+1:N] = [(x[1]-particles[i][1])/(x[0]-particles[i][0]) for x in particles[i+1:]]
    mat_dist = mat_dist + mat_dist.transpose()
    mat_dist=np.exp(mat_dist*2j)
    
    mat = mat_hank_zero + mat_dist*mat_hank_two
    #Correct the value for the main diagonal
    np.fill_diagonal(mat,1)
       
    
    return mat*(-Gamma_1/2)


def light_perp_2(mat,particles,Gamma_1,N):
    #Create matrix with coordinates distance and sum the components
    mat_dist = zeros((N,N))
    #Here phi_jl = (y[j]-y[l])/(x[j]-x[l]), and so, phi_jl=phi_lj
    #Also, the line 'i' must contain phi_i0,phi_i1, ..., phi_iN
    #So, on line 'i' the particle 'i' coordinates remain fixed, while the others vary for all particles in the system
    
    for i in range(0, N):
        mat_dist[i, i+1:N] = [(x[1]-particles[i][1]) / (x[0]-particles[i][0]) if x[0] != particles[i][0] else 0 for x in particles[i+1:]]

    mat_dist = mat_dist + mat_dist.transpose()
    mat_dist=np.exp(mat_dist*2j)
    print(np.max(mat_dist))
    #Perform Hankel function for all the matrix
    #Add 1 to the main diagonal just to don't crash the vectorized calculus
    fill_diagonal(mat, 1)
    mat_hank_zero=H0(mat)
    mat_hank_two = H2(mat)
    print(np.max(mat_hank_zero))
    print(np.max(mat_hank_two))
    
    mat = mat_hank_zero + mat_dist*mat_hank_two
    print(np.max(mat_dist))
    #Correct the value for the main diagonal
    np.fill_diagonal(mat,1)
       
    
    return mat*(-Gamma_1/2)