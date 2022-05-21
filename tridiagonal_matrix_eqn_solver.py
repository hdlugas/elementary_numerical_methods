# this script solves Ax=b for tridiagonal matrix A
#exercise 2.3.C1

import numpy as np


def tri(n,a,d,c,b):
    #a,d,c are nx1 arrays and subdiagonal, main diagonal, and superdiagonal, respectively
    x = np.zeros(n)
    count = 0
    for i in range(1,n):
        xmult = a[i-1]/d[i-1]
        d[i] = d[i]-xmult*c[i-1]
        b[i] = b[i]-xmult*b[i-1]
        count += 3
    x[n-1] = b[n-1]/d[n-1]
    count += 1
    for i in range(n-2,-1,-1):
        x[i] = (b[i]-c[i]*x[i+1])/d[i]
        count += 2
    return(x,count)


n = 3
a = [2,6,0]
d = [1,4,7]
c = [3,5,0]
b = [0,1,2]
x,count = tri(n,a,d,c,b)

A=np.asarray([[1,3,0],[2,4,5],[0,6,7]])
x_np = np.linalg.solve(A,b)


print('tridiagonal matrix solver solution: x =',x)
print('numpy\'s solution:                 : x =',x_np)
print('the number of multiplication operations used to compute x is ',count)
