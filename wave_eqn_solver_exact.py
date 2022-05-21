#this script approximates the solution of the wave equation PDE using an odd periodic extension of the solution function

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def f(x):
    return x**2


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

#test the odd periodic extension function
h = pow(2,-5)
k = pow(2,-6)
n = int((1/h))
m = int((1/k))


domain = np.linspace(-5,5,100) #input only odd numbers for the bounds, otherwise have to condition on y%2==0,1 in the for loop below
y = np.zeros((len(domain)))
for i in range(0,len(y)):
    if i == 0 or i == len(y)-1:
        y[i]=f(1)
    else:
        y[i] = eval_odd_periodic_ext(f,domain[i],h)


plt.scatter(domain,y,marker='.')
plt.xlabel('x - axis')
plt.ylabel('y - axis')
plt.title('Odd Periodic Extension of f(x) with original domain [0,1]')
plt.show()




#now test the wave equation solver
solution = wave_eqn_solver_exact(h,k,n,m,f,eval_odd_periodic_ext)

x = np.linspace(0,1,n)
y = np.linspace(0,7,m)
x,y = np.meshgrid(x,y)
fig,ax = plt.subplots(subplot_kw={"projection": "3d"})
surf_exact = ax.plot_surface(x,y,solution,linewidth=0,antialiased=False,color='blue')
ax.set_title('Solution Surface')
ax.set_xlabel('Space')
ax.set_ylabel('Time')
plt.show()







