#this script solves a pentadiagonal matrix

import numpy as np


def penta(n,e,a,d,c,f,b):
    #from bottom to top, diagonals are e,a,d,c,f, all nx1 arrays
    #the last element of a and c and the last two elements of e and f doesn't matter
    x = np.zeros(n)
    count = 0
    r = a[0]
    s = a[1]
    t = e[0]
    for i in range(1,n-1):
        xmult = r/d[i-1]
        d[i] -= xmult*c[i-1]
        c[i] -= xmult*f[i-1]
        b[i] = b[i]-xmult*b[i-1]
        xmult = t/d[i-1]
        r = s-xmult*c[i-1]
        d[i+1] -= xmult*f[i-1]
        b[i+1] = b[i+1]-xmult*b[i-1]
        s = a[i+1]
        t = e[i]
        count = count+8
    xmult = r/d[n-2]
    d[n-1] -= xmult*c[n-2]
    x[n-1] = (b[n-1]-xmult*b[n-2])/d[n-1]
    x[n-2] = (b[n-2]-c[n-2]*x[n-1])/d[n-2]
    count += 6
    for i in range(n-3,-1,-1):
        x[i] = (b[i]-f[i]*x[i+2]-c[i]*x[i+1])/d[i]
        count += 3
    return (x,count)


n = 5
e = [1,2,1,0,0]
a = [1,1,-1,1,0]
d = [2,-2,1,2,1]
c = [-1,1,2,1,0]
f = [1,1,-1,0,0]
b = [1,2,1,4,4]

A=[[2,-1,1,0,0],[1,-2,1,1,0],[1,1,1,2,-1],[0,2,-1,2,1],[0,0,1,1,1]]

x,count = penta(n,e,a,d,c,f,b)
print('\nThe solution using the pentadiagonal solver is x =',x)
print('\nThe solution using numpy\'s matrix solver is x =',np.linalg.solve(A,b))
print('\nThe number of long arithmetic operations used was ',count,'\n')
