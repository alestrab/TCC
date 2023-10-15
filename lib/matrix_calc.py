from numpy import fill_diagonal,zeros,ones
from hankel_func import *

def light_par(mat,Gamma_0):
    #Perform Hankel function for all the matrix
    #Add 1 to the main diagonal just to don't crash the vectorized calculus
    fill_diagonal(mat, 1)
    mat = H0(mat)

    #Correct the value for the main diagonal
    fill_diagonal(mat, 1)
    mat = mat*(-Gamma_0/2)
    return mat

def light_perp(mat,particles,Gamma_1,N):
    #Perform Hankel function for all the matrix
    #Add 1 to the main diagonal just to don't crash the vectorized calculus
    fill_diagonal(mat, 1)
    mat = H0(mat)
    
    mat_dist = zeros((N,N))
    for i in range(0,N):
        mat_dist[i,i+1:N] = [(x[1]-particles[i][1])/(x[0]-particles[i][0]) for x in particles[i+1:]]
    mat_dist = mat_dist + mat_dist.transpose()+ones((N,N))
       
    #Correct the value for the main diagonal
    fill_diagonal(mat, 1)
    mat = mat*mat_dist*(-Gamma_1/2)
    
    return mat