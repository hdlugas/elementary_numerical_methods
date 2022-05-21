

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def tri(n,a,d,c,b):
    #a,d,c are nx1 arrays and subdiagonal, main diagonal, and superdiagonal, respectively
    x = np.zeros(n)
    for i in range(1,n):
        xmult = a[i-1]/d[i-1]
        d[i] = d[i]-xmult*c[i-1]
        b[i] = b[i]-xmult*b[i-1]
    x[n-1] = b[n-1]/d[n-1]
    for i in range(n-2,-1,-1):
        x[i] = (b[i]-c[i]*x[i+1])/d[i]
    return(x)



def parabolic_pde_solver_crank_nicolson(h,k,n,m):
    s = h*h/k
    r = s+2
    d = [r]*(n-1)
    c = [-1]*(n-1)
    u = np.zeros((n-1))
    v = np.zeros((n-1))
    w = np.zeros((n-1))
    sol_approx = np.zeros((m-1,n-1))#top row is solution at t_0, then time step increases as you move down columns
    sol_exact = np.zeros((m-1,n-1))
    for i in range(1,n):
        u[i-1] = np.sin(np.pi*i*h)
    sol_approx[0,:] = u
    sol_exact[0,:] = u
    for j in range(1,m-1):
        d = [r]*(n-1)
        for i in range(0,n-1):
            v[i] = u[i]*s
        v = tri(n-1,c,d,c,v)
        sol_approx[j,:] = v
        t = j*k
        for i in range(1,n):
            w[i-1] = pow(np.e,-1*np.pi*np.pi*t)*np.sin(np.pi*i*h)
        sol_exact[j,:] = w
        u=v
    return sol_approx,sol_exact


h = pow(2,-4)
k = pow(2,-10)

n = int((1/h))

m = int((0.25/k))


sol_approx,sol_exact = parabolic_pde_solver_crank_nicolson(h,k,n,m)
x = np.linspace(0,1,n-1)
y = np.linspace(0,0.25,m-1)
x,y = np.meshgrid(x,y)


fig,ax = plt.subplots(subplot_kw={"projection": "3d"})
surf_approx = ax.plot_surface(x,y,sol_approx,linewidth=0,antialiased=False,color='blue')
ax.set_title('Approximate Solution Surface')
ax.set_xlabel('Space')
ax.set_ylabel('Time')
plt.show()


fig,ax = plt.subplots(subplot_kw={"projection": "3d"})
surf_exact = ax.plot_surface(x,y,sol_exact,linewidth=0,antialiased=False,color='blue')
ax.set_title('Exact Solution Surface')
ax.set_xlabel('Space')
ax.set_ylabel('Time')
plt.show()


fig,ax = plt.subplots(1,1)
ax.contour(x,y,sol_approx,levels=50)
ax.contour(x,y,sol_exact,levels=50)
ax.set_title('Contour Plot of Exact vs Approximate Solution')
ax.set_xlabel('Space')
ax.set_ylabel('Time')
plt.show()
