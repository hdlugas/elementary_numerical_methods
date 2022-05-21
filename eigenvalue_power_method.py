#this script uses the power method to calculate the dominant eigenvalue of a matrix
#exercise 8.3.C1

import numpy as np
from scipy import linalg as la


def phi(x):
    return x[0]


def power_method(A,phi,x,n_iter):
    y = np.zeros(n_iter)
    for i in range(0,n_iter):
        y = np.matmul(A,x)
        r = phi(y)/phi(x) #r approaches the dominant eigenvalue as n_iter->infinity
        x = y
    return r


A = np.zeros((4,4))
A = [[5,4,1,1],[4,5,1,1],[1,1,4,2],[1,1,2,4]]
print('\n Test matrix:\n',A)

eig_val = la.eigvals(A)
mag_vec = np.zeros(len(eig_val))
for i in range(0,len(eig_val)):
    mag_vec[i] = np.abs(eig_val[i])
print('\nUsing scipy\'s function to calculate the eigenvalues of this matrix, we obtain:\n', eig_val)

n_iter = 10
x = np.asarray([1,2,3,4])
eig = power_method(A,phi,x,n_iter)
print('\nOur computed dominant eigenvalue using ',n_iter,' iterations of the power method is ',eig)
print('\nThe absolute error between our computed eigenvalue and scipy\'s computed eigenvalue is ',np.abs(eig-max(eig_val)),'\n')
