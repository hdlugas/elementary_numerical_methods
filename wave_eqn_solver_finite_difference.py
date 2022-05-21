#this script approximates the solution to the wave equation PDE using two different finite difference methods

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.linalg as la

def f(x):
    return 0.5-np.abs(x-0.5)

def eval_odd_periodic_ext(f,x,h): #evaluate g(x) where g is the odd periodic extension of f with f's domain being [0,1]
    iter_x = int((1/h))
    g = np.linspace(-5,5,10*iter_x)
    if x >= 0 and x <= 1:
        g = f(x)
    if x < 0 and x >= -1:
        g = -f(-x)
    if x > 1 and np.floor(x)%2 == 0:
        g = f(x-np.floor(x))
    if x > 1 and np.floor(x)%2 == 1:
        g = -f(np.ceil(x)-x)
    if x < -1 and np.ceil(x)%2 == 0:
        g = -f(np.ceil(x)-x)
    if x < -1 and np.ceil(x)%2 == 1:
        g = f(x-np.floor(x))
    return g


def wave_eqn_solver_exact(h,k,n,m,f,eval_odd_periodic_ext):
    u = np.zeros((n))
    solution = np.zeros((m,n))#top row is solution at t_0, then time step increases as you move down columns
    for j in range(0,m):
        t = j*k
        for i in range(0,n):
            x = i*h
            u[i] = 0.5*(eval_odd_periodic_ext(f,x+t,h)+eval_odd_periodic_ext(f,x-t,h))
        solution[j,:] = u
    return solution


def hyperbolic_pde_solver1(h,k,n,m,f,eval_odd_periodic_ext):
    u = np.zeros((n+1))
    v = np.zeros((n+1))
    w = np.zeros((n+1))
    rho = pow(k/h,2)
    solution = np.zeros((m,n+1))#top row is solution at t_0, then time step increases moving down columns
    for i in range(1,n):
        x = i*h
        w[i] = eval_odd_periodic_ext(f,x,h)
        v[i] = eval_odd_periodic_ext(f,x,h)
    solution[0,:] = v
    for j in range(1,m):
        for i in range(1,n):
            u[i] = rho*(v[i+1]+v[i-1])+2*(1-rho)*v[i]-w[i]
        solution[j,:] = u
        for i in range(1,n):
            w[i] = v[i]
            v[i] = u[i]
    return solution



def hyperbolic_pde_solver2(h,k,n,m,f,eval_odd_periodic_ext):
    u = np.zeros((n+1))
    v = np.zeros((n+1))
    w = np.zeros((n+1))
    rho = pow(k/h,2)
    solution = np.zeros((m,n+1))#top row is solution at t_0, then time step increases moving down columns
    for i in range(1,n):
        x = i*h
        w[i] = eval_odd_periodic_ext(f,x,h)
        v[i] = 0.5*rho*(eval_odd_periodic_ext(f,x-h,h)+eval_odd_periodic_ext(f,x+h,h))+(1-rho)*eval_odd_periodic_ext(f,x,h)
    solution[0,:] = v
    for j in range(1,m):
        for i in range(1,n):
            u[i] = rho*(v[i+1]+v[i-1])+2*(1-rho)*v[i]-w[i]
        solution[j,:] = u
        for i in range(1,n):
            w[i] = v[i]
            v[i] = u[i]
    return solution


h = 1/16
k = 1/32
n = int((1/h))
m = int((0.4/k))

solution1 = hyperbolic_pde_solver1(h,k,n,m,f,eval_odd_periodic_ext)
solution2 = hyperbolic_pde_solver2(h,k,n,m,f,eval_odd_periodic_ext)
solution_exact = wave_eqn_solver_exact(h,k,n+1,m,f,eval_odd_periodic_ext)

x = np.linspace(0,1,n+1)
y = np.linspace(0,0.4,m)
x,y = np.meshgrid(x,y)
x_exact = np.linspace(0,1,n+1)
y_exact = np.linspace(0,0.4,m)
x_exact,y_exact = np.meshgrid(x_exact,y_exact)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax = fig.gca(projection='3d')
ax.set_title('Solution Surfaces of Exact (orange) and Approximate (blue) Solutions using method 1')
ax.set_xlabel('Space')
ax.set_ylabel('Time')
surf1 = ax.plot_surface(x,y,solution1,color='blue')
surf_exact = ax.plot_surface(x_exact,y_exact,solution_exact,color='orange')

fig = plt.figure()
ax = plt.axes(projection='3d')
ax = fig.gca(projection='3d')
ax.set_title('Solution Surfaces of Exact (orange) and Approximate (blue) Solutions using method 2')
ax.set_xlabel('Space')
ax.set_ylabel('Time')
surf2 = ax.plot_surface(x,y,solution2,color='blue')
surf_exact = ax.plot_surface(x_exact,y_exact,solution_exact,color='orange')


fig,ax = plt.subplots(1,1)
ax.contour(x,y,solution1,levels=50)
ax.contour(x_exact,y_exact,solution_exact,levels=50)
ax.set_title('Contour Plot of Exact vs Approximate solutions using method 1')
ax.set_xlabel('Space')
ax.set_ylabel('Time')
plt.show()

fig,ax = plt.subplots(1,1)
ax.contour(x,y,solution2,levels=50)
ax.contour(x_exact,y_exact,solution_exact,levels=50)
ax.set_title('Contour Plot of Exact vs Approximate solutions using method 2')
ax.set_xlabel('Space')
ax.set_ylabel('Time')
plt.show()


norm1 = la.norm(solution1-solution_exact,ord='fro')
norm2 = la.norm(solution2-solution_exact,ord='fro')
print('\nThe Frobenius norm of the difference between the exact and our approximate solution using method 1 is',norm1)
print('\nThe Frobenius norm of the difference between the exact and our approximate solution using method 2 is',norm2)

