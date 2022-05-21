#this script tests the following integration rules:
#composite trapezoidal rule
#composite Simpson's rule
#two-point Gaussian quadrature rule
#three-point Gaussian quadrature rule

import numpy as np
import matplotlib.pyplot as plt




def f(x):
    return pow(np.e,-1*x*x)


def trap(a,b,f,I):
    print('\nTrapezoid Method:')
    n = 1
    error = []
    for i in range(1,8):
        h = (b-a)/n
        sum = 0.5*(f(a)+f(b))
        for i in range(1,n):
            x = a+h*i
            sum += f(x)
        sum = sum*h
        print('\nThe given value of the integral is                   ',I)
        print('Using ',n,'subintervals, the approximate value of the integral is',sum)
        print('The absolute error is ',np.abs(I-sum),'\n')
        error.append(np.abs(I-sum))
        n = n*2
    return error





def simp(a,b,f,I):
    print('\n\n\n\nSimpson Method:')
    n = 1
    error = []
    for i in range(1,8):
        h = (b-a)/n
        sum = f(a)+f(b)
        sum1 = 0
        sum2 = 0
        if n > 1:
            for i in range(1,int(n/2)+1):
                x = a+(2*i-1)*h
                sum1 += f(x)
            if n > 3:
                for i in range(1,int((n-2)/2)+1):
                    x = a+2*i*h
                    sum2 += f(x)
            sum += 4*sum1+2*sum2
            sum += h/3
        print('\nThe given value of the integral is             ',I)
        print('Using ',n,'subintervals, the approximate value of the integral is',sum)
        print('The absolute error is ',np.abs(I-sum),'\n')
        error.append(np.abs(I-sum))
        n = n*2
    return error


def lamb(a,b,x):
    return ((b-a)/2)*x+a/2+b/2



def gq2p(a,b,f,I,lamb):
    print('\n\n\n\nTwo-point Gaussian Quadrature:')
    n = 1
    a0 = a
    b0 = b
    error = []
    for i in range(1,8):
        sum = 0
        h = (b0-a0)/n
        endpoints = np.linspace(a0,b0,n+1)
        for j in range(0,len(endpoints)-1):
            a = endpoints[j]
            b = endpoints[j+1]
            sum += ((b-a)/2)*(5/9*f(lamb(a,b,-1*np.sqrt(3/5)))+8/9*f(lamb(a,b,np.sqrt(3/5))))

        print('\nThe given value of the integral is             ',I)
        print('Using ',n,'subintervals, the approximate value of the integral is',sum)
        print('The absolute error is ',np.abs(I-sum),'\n')
        error.append(np.abs(I-sum))
        n = n*2
    return error




def gq4p(a,b,f,I,lamb):
    print('\n\n\n\nFour-point Gaussian Quadrature:')
    n = 1
    a0 = a
    b0 = b
    error = []
    for i in range(1,8):
        sum = 0
        h = (b0-a0)/n
        endpoints = np.linspace(a0,b0,n+1)
        for j in range(0,len(endpoints)-1):
            a = endpoints[j]
            b = endpoints[j+1]
            sum += ((b-a)/2)*((1/2+np.sqrt(10/3)/12)*f(lamb(a,b,-1*np.sqrt(1/7+(3-4*np.sqrt(0.3)))))+
            (1/2-np.sqrt(10/3)/12)*f(lamb(a,b,-1*np.sqrt(1/7*(3+4*np.sqrt(0.3)))))+
            (1/2+np.sqrt(10/3)/12)*f(lamb(a,b,np.sqrt(1/7*(3-4*np.sqrt(0.3)))))+
            (1/2-np.sqrt(10/3)/12)*f(lamb(a,b,np.sqrt(1/7*(3+4*np.sqrt(0.3))))))

        print('\nThe given value of the integral is             ',I)
        print('Using ',n,'subintervals, the approximate value of the integral is',sum)
        print('The absolute error is ',np.abs(I-sum),'\n')
        error.append(np.abs(I-sum))
        n = n*2
    return error




a = 0
b = 1
I = 0.7468241328124270

n_vec = []
n = 1
for i in range(1,8):
    n_vec.append(n)
    n = n*2
x = np.log2(n_vec)



error_t = trap(a,b,f,I)
y1 = np.log10(error_t)
t = 2
for i in range(1,len(error_t)):
    print('error __',t,'=',error_t[i],'    log_2(e_{n/2}/e_{n}) =',np.log2(error_t[i-1]/error_t[i]))
    t = t*2


error_s = simp(a,b,f,I)
y2 = np.log10(error_s)
t = 2
for i in range(1,len(error_s)):
    print('error __',t,'=',error_s[i],'    log_2(e_{n/2}/e_{n}) =',np.log2(error_s[i-1]/error_s[i]))
    t = t*2



error_gq2p = gq2p(a,b,f,I,lamb)
y3 = np.log10(error_gq2p)
t = 2
for i in range(1,len(error_gq2p)):
    print('error __',t,'=',error_gq2p[i],'    log_2(e_{n/2}/e_{n}) =',np.log2(error_gq2p[i-1]/error_gq2p[i]))
    t = t*2



error_gq4p = gq4p(a,b,f,I,lamb)
y4 = np.log10(error_gq4p)
t = 2
for i in range(1,len(error_gq4p)):
    print('error __',t,'=',error_gq4p[i],'    log_2(e_{n/2}/e_{n}) =',np.log2(error_gq4p[i-1]/error_gq4p[i]))
    t = t*2


plt.plot(x,y1)
plt.plot(x,y2)
plt.plot(x,y3)
plt.plot(x,y4)
plt.xlabel('log_2(n)')
plt.ylabel('log_10(e_{n/2}/e_{n}')
plt.title('log-log Plot')
plt.legend(["Trapezoid Method","Simpson\'s Method", "Two-Point Gaussian Quadrature Rule","Four-Point Gaussian Quadrature Rule"],loc ="upper right")
plt.show()









