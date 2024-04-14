import math
import numpy as np
from matplotlib import pyplot as plt

def K (x, t):
    return 1 / 2 * x * t

def f (x):
    return 5 / 6 * x

def u (x):
    return x   

def u_n (x, h, U):
    u = f(x)
    for i in range (len(U)):
        u += h * K(x, i * h) * U[i]
    return u

def solve (h, n):

    x_k = [0] * (n - 1)
    A = np.zeros((n - 1, n - 1))
    B = np.zeros((n - 1, 1))

    for i in range (n - 1):
        x_k[i] = h * (i + 1)
        
    for i in range (n - 1):
        A[i][i] = 1    
        B[i][0] = f(x_k[i])
        for j in range (n - 1):
            A[i][j] += -1 * h * K(x_k[i], x_k[j])
    U = np.linalg.solve(A, B)
#    print(U)
    return U

def dist1 (x, h1, U1, h2, U2):
    return pow((u_n(x, h1, U1) - u_n (x, h2, U2)), 2)

def dist2 (x, h1, U1, h2, U2):
    return pow((u(x) - u_n (x, h2, U2)), 2)

def Gauss3 (a, b, func, h1, U1, h2, U2):        
        x = [(b + a) / 2 - 0.7745967 * (b - a) / 2, (b + a) / 2, (b + a) / 2 + 0.7745967 * (b - a) / 2]
        c = [0.5555556, 0.8888889, 0.5555556]
        
        t = 0
        for i in range (3):
            t += c[i] * func(x[i], h1, U1, h2, U2)
        return t * (b - a) / 2
    
def Gauss5 (a, b, func, h1, U1, h2, U2):
        x = [(b + a) / 2 - 0.9061798 * (b - a) / 2, (b + a) / 2 - 0.5384693 * (b - a) / 2,(b + a) / 2, (b + a) / 2 + 0.5384693 * (b - a) / 2,(b + a) / 2 + 0.9061798 * (b - a) / 2]
        c = [0.2369269, 0.4786287, 0.5688888, 0.4786287, 0.2369269]
        
        t = 0
        for i in range (5):
            t += + c[i] * func(x[i], h1, U1, h2, U2)
        return t * (b - a) / 2 

def algorithm (U, U1, h, func):
    alpha = a
    beta = b
    I = 0
    E = 0
    k = 0
    k_max = 0
    while (alpha < beta):
        h1 = beta - alpha
        I_1 = Gauss3(alpha, beta, func, 2 * h, U, h, U1)
        I_2 = Gauss5(alpha, beta, func, 2 * h, U, h, U1)
        delta = I_1 - I_2
        if abs(delta) < max(eps_A, eps_O * abs(I)):
            k = 0
            alpha = beta
            beta = b
            I = I + I_2 + delta
            E = E + delta
        else:
            beta = (alpha + beta) / 2
            k = k + 1
            if k > k_max:
                k = 0
#                print("feature point: ", alpha, beta)
                I = I_2 + delta
                E = E + delta
                alpha = beta
                beta = b
#    print (I, E)
    return math.sqrt(I)
        
    
n = 64
eps_A = 0.000001
eps_O = 0.001
a = 0
b = 1
h = (b - a) / n

U = solve(h, n)
norm = 1
i = 0
while (norm > eps_O):
    i += 1
    n *= 2
    h = (b - a) / n
    U1 = solve(h, n)
    norm = algorithm (U, U1, h, dist1)
    U = U1
    print ("n: ", n, "norm: ", norm)
print (i, " step")
print ("|u - u_n| = ", algorithm(U, U1, h, dist2))
    
    
#print graph

#plt.plot(x_k, U)
#plt.show()


    
    
