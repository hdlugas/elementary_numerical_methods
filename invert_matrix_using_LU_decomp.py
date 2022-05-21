#this script inverts a nonsingular matrix using LU decomposition

import numpy as np
from scipy import linalg as la

def invert_matrix(L,U):
    n = len(L[:,0])
    x = []
    y = []
    for i in range(0,n):
        b = np.zeros(n)
        b[i] = 1
        y = la.solve(L,b)
        sln_vec = la.solve(U,y)
        x.append(sln_vec)

    A = np.zeros((n,n))
    for i in range(0,n):
        for j in range(0,n):
            A[i,j] = x[j][i]
    return(A)


def LU_decomp(LU):
    n = len(LU[:,0])
    L = np.zeros((n,n))
    U = np.zeros((n,n))
    for i in range(0,n):
        for j in range(0,n):
            if i == j:
                L[i,j] = 1
            if i > j:
                L[i,j] = LU[i,j]
            if i <= j:
                U[i,j] = LU[i,j]
    return(L,U)




n = 3
A = np.zeros((n,n))
L = np.zeros((n,n))
U = np.zeros((n,n))
A[0][0] = 2
A[0][1] = 3 
A[0][2] = 0
A[1][0] = 4
A[1][1] = 5
A[1][2] = 8
A[2][0] = 0 
A[2][1] = 6
A[2][2] = 7

print('\nTest matrix:\n',A)
LU,piv = la.lu_factor(A)
L,U = LU_decomp(LU)

sln = invert_matrix(L,U)
print('\nA^{-1} can be found by rearranging the rows and columns of the matrix:\n',sln,'\n')

#print('Check our solution:\n')
#print('A^{-1}A =',np.matmul(sln,A))
#print('AA^{-1} =',np.matmul(A,sln))
