#this script uses the Jacobi, Gauss-Seidel, and Successive Over-Relaxatoin methods to iteratively solve a linear system of equations


import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la



def Jacobi(A,b,x,epsilon,sln_vec):
    n = len(b)
    n_jacobi = 0
    residual_norm_vec = []
    while np.abs(la.norm(sln_vec,ord=2)-la.norm(x,ord=2)) > epsilon:
        n_jacobi += 1
        y = x
        x = np.zeros((n))
        for i in range(0,n):
            sum = b[i]
            for j in range(0,n):
                if i != j:
                    sum = sum-A[i,j]*y[j]
            x[i] = sum/A[i,i]
        residual_norm_vec.append(la.norm(x,ord=2))
    return x,n_jacobi,residual_norm_vec


def gauss_seidel(A,b,x,epsilon,sln_vec):
    n = len(b)
    n_gs = 0
    residual_norm_vec = []
    while np.abs(la.norm(sln_vec,ord=2)-la.norm(x,ord=2)) > epsilon:
        n_gs += 1
        for i in range(0,n):
            sum = b[i]
            for j in range(0,n):
                if i != j:
                    sum = sum-A[i,j]*x[j]
            x[i] = sum/A[i,i]
        residual_norm_vec.append(la.norm(x,ord=2))
    return x,n_gs,residual_norm_vec



def sor(A,b,x,w,epsilon,sln_vec):
    n = len(b)
    n_sor = 0
    residual_norm_vec = []
    while np.abs(la.norm(sln_vec,ord=2)-la.norm(x,ord=2)) > epsilon:
        n_sor += 1
        y = x
        for i in range(0,n):
            sum = b[i]
            for j in range(0,n):
                if i != j:
                    sum = sum-A[i,j]*x[j]
            x[i] = sum/A[i,i]
            x[i] = w*x[i]+(1-w)*y[i]
        residual_norm_vec.append(la.norm(x,ord=2))
    return x,n_sor,residual_norm_vec


n = 25
d = [2]*n
c = [-1]*n
a = [-1]*n
b = [1]*n
x0 = [0]*n
epsilon = 0.0005
w = 1.1

A=np.zeros((n,n))
for i in range(0,n):
    for j in range(0,n):
        if i == j+1:
            A[i,j] = -1
        if i == j-1:
            A[i,j] = -1
        if i == j:
            A[i,j] = -2
print('\nTest matrix for n =',n,'case:\n',A)


sln_vec = la.solve(A,b)
sln_j,n_jacobi,residual_norm_vec_jacobi = Jacobi(A,b,x0,epsilon,sln_vec)
x0 = [0]*n
sln_gs,n_gs,residual_norm_vec_gs = gauss_seidel(A,b,x0,epsilon,sln_vec)
x0 = [0]*n
sln_sor,n_sor,residual_norm_vec_sor = sor(A,b,x0,w,epsilon,sln_vec)


print('\nSolution using the Jacobi method:\n',sln_j)
print('\nSolution using the Gauss-Seidel method:\n',sln_gs)
print('\nSolution using the Successive Overrelaxation method:\n',sln_sor)
print('\nScipy\'s solution:\n',sln_vec)

print('\nThe number of iterations needed for the euclidean norm of the solution vector from the Jacobi')
print('method to be within +- ',epsilon,'of the euclidean norm of the solution vector is ',n_jacobi)

print('\nThe number of iterations needed for the euclidean norm of the solution vector from the Gauss-Seidel')
print('method to be within +- ',epsilon,'of the euclidean norm of the solution vector is ',n_gs)

print('\nThe number of iterations needed for the euclidean norm of the solution vector from the Successive Overrelaxation')
print('method to be within +- ',epsilon,'of the euclidean norm of the solution vector is ',n_sor,'\n')

t_j = np.zeros(len(residual_norm_vec_jacobi))
for i in range(0,len(residual_norm_vec_jacobi)):
    t_j[i] = i+1

t_gs = np.zeros(len(residual_norm_vec_gs))
for i in range(0,len(residual_norm_vec_gs)):
    t_gs[i] = i+1

t_sor = np.zeros(len(residual_norm_vec_sor))
for i in range(0,len(residual_norm_vec_sor)):
    t_sor[i] = i+1

plt.scatter(t_j,np.log10(residual_norm_vec_jacobi),marker='.',label='Jacobi Method')
plt.scatter(t_gs,np.log10(residual_norm_vec_gs),marker='.',label='Gauss-Seidel Method')
plt.xlabel('Iteration number')
plt.ylabel('Logarithm base 10 of euclidean norm of residual vector')
plt.title('Convergence Comparison of Jacobi, Gauss-Seidel,\nand Successive Overrelaxation Methods')
plt.legend()
plt.show()

plt.scatter(t_j,np.log10(residual_norm_vec_jacobi),marker='.',label='Jacobi Method')
plt.scatter(t_sor,np.log10(residual_norm_vec_sor),marker='.',label='Successive Overrelaxation Method')
plt.xlabel('Iteration number')
plt.ylabel('Logarithm base 10 of euclidean norm of residual vector')
plt.title('Convergence Comparison of Jacobi, Gauss-Seidel,\nand Successive Overrelaxation Methods')
plt.legend()
plt.show()

plt.scatter(t_gs,np.log10(residual_norm_vec_gs),marker='.',label='Gauss-Seidel Method')
plt.scatter(t_sor,np.log10(residual_norm_vec_sor),marker='.',label='Successive Overrelaxation Method')
plt.xlabel('Iteration number')
plt.ylabel('Logarithm base 10 of euclidean norm of residual vector')
plt.title('Convergence Comparison of Jacobi, Gauss-Seidel,\nand Successive Overrelaxation Methods')
plt.legend()
plt.show()



print('\nSince the vector of residual norms of the Gauss-Seidel and SOR methods seem to be nearly identical from the graph,')
print('consider their distance: ',np.asarray(residual_norm_vec_gs)-np.asarray(residual_norm_vec_sor))

