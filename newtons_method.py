#this script uses Newton's method to compute the reciprocal of a number R using
#the iterative formula x_{n+1}=x_{n}(2-x_{n}R)
#exercise 3.2.23


import numpy as np
import matplotlib.pyplot as plt

def Newton(f,df,R,x,nmax,epsilon,delta):
    fx = f(x,R)
    for i in range(0,nmax-1):
        fp = df(x,R)
        if np.abs(fp) < delta:
            print('small derivative, try new starting point')
        x = x-fx/fp
        fx = f(x,R)
        print('iteration: ',i+1)
        print('x =',x)
        print('R =',R)
        print('|x-1/R|   =',np.abs(x-1/R))
        print('|x-1/R|^2 =',np.abs(x-1/R)**2)
        print('epsilon =',epsilon)
        print('f(x) =',fx,'\n')
        if np.abs(x-1/R) < epsilon:
            print('Since |x-1/R|<epsilon, the algorithm has converged at x =',x,'.\n')
            print('Also notice the quadratic convergence exhibited, i.e. |x_{n+1}-1/R| ~ c*|x_{n}-1/R|^2 for c=1')
            print('We notice that this is the case since |x_{n+1}-1/R| and |x_{n}-1/R|^2 are about the same order of magnitude for each n.\n')
            break


def f(x,R):
        y = R-1/x
        return y

def df(x,R):
        y = 1/(pow(x,2))
        return y


nmax = 600
x0 = 0.2
R = 4
epsilon = 5*10**(-7)
delta = 10**(-8)
Newton(f,df,R,x0,nmax,epsilon,delta)


x = np.linspace(0.1,1,100)
f = np.zeros(100)
g = np.zeros(100)

for i in range(0,len(x)):
    f[i] = 4-1/x[i]
    g[i] = 0


plt.scatter(x,f,label='f(x)=4-1/x')
plt.scatter(x,g,label='g(x)=0')
plt.xlabel('x - axis')
plt.ylabel('y - axis')
plt.legend()
plt.show()
