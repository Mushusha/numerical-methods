import math
import numpy as np

def func (x):
    return math.exp(pow((x[0] - 1), 2) + 2 * pow((x[1] - 2), 2))

def grad (x):
    gradient = [0] * len(x)
    for i in range (len(x)):
        x[i] = x[i] + h
        gradient[i] = func(x)
        x[i] = x[i] - 2 * h
        gradient[i] = gradient[i] - func(x)
        x[i] = x[i] + h
        gradient[i] = gradient[i] / (2 * h)
    return gradient

def tFunc (x, t):
    temp = [0] * n
    for i in range (n):
        temp[i] = x[i] - t * MNeuton(x)[i]
        
    return func(temp)

def tGrad (x, t):
    gradient = tFunc(x, t + h)
    gradient = gradient - tFunc(x, t - h)
    gradient = gradient / (2 * h)

    return gradient


def cond1 (x1, x2, eps):
    norm = 0
    for i in range (len(x1)):
        norm = norm + pow(x2[i] - x1[i], 2)
    if (math.sqrt(norm) <= eps):
        return 1
    return 0

def cond2 (x1, x2, eps):
    if (abs(func(x2) - func(x1)) <= eps):
        return 1
    return 0

def cond3 (x, eps):
    norm = 0
    for i in range (len(x)):
        norm = norm + pow(grad(x)[i], 2)
    if (math.sqrt(norm) <= eps):
        return 1
    return 0

def cond (x1, x2, eps):
    if (cond1(x1, x2, eps) & cond2(x1, x2, eps) & cond3(x2, eps)):
        return 1
    return 0

def c (x, a, b):
    return (tFunc(x, a) - tFunc(x, b) + tGrad(x, b) * b - tGrad(x, a) * a) / (tGrad(x, b) - tGrad(x, a))
    

def MNeuton (x):
    n = len(x)
    d = np.zeros((n, n))
    g = np.zeros((n, 1))
    ans = [0] * n
    for i in range (n):
        for j in range (n):
            if i == j:
                x[i] = x[i] - h
                d[i][j] = func(x)
                x[i] = x[i] + h
                d[i][j] = d[i][j] - 2 * func(x)
                x[i] = x[i] + h
                d[i][j] = (d[i][j] + func(x)) / (h * h)
                x[i] = x[i] - h
            if i != j:
                d[i][j] = d[i][j] - func(x)
                x[i] = x[i] + h
                d[i][j] = d[i][j] + func(x)
                x[j] = x[j] - h
                d[i][j] = d[i][j] - func(x)
                x[i] = x[i] - h
                d[i][j] = (d[i][j] + func(x)) / (h * h)
                x[j] = x[j] + h
        g[i] = grad(x)[i]
        
    return (np.matmul(np.linalg.inv(d), g))
                



n = 2
h = 0.01
x0 = [0] * n
eps = 0.00001
x0[0] = 1.1
x0[1] = 1.9

 # 1 method

k = 1
x1 = [0] * n
x2 = [0] * n
beta = 3
l = 0.5
m = 2
alpha = 1
for i in range(n):
    x2[i] = x0[i]
    x1[i] = x0[i]

#    x2[i] = x1[i] - a * grad(x1, a)[i]
while (cond(x1, x2, math.sqrt(eps)) != 1):
    alpha = beta
    y = [0] * n
    z = [0] * n
    for i in range (n):
                y[i] = x2[i] - alpha * grad(x2)[i]
                z[i] = x2[i] - beta * grad(x2)[i]    
    if (alpha == beta):
        while func(y) >= func(x2):
            alpha = l * alpha
            for i in range (n):
                y[i] = x2[i] - alpha * grad(x2)[i]

            
    if (alpha == beta):
        while func(y) < func(z):
            alpha = m * alpha
            for i in range (n):
                y[i] = x2[i] - alpha * grad(x2)[i]
                z[i] = x2[i] - beta * grad(x2)[i]

    for i in range(n):
        x1[i] = x2[i]
        x2[i] = x1[i] - alpha * grad(x1)[i]
    k += 1
    
    
print(k, x2)

 # 2 method
 
a = 0
b = 1
k = 0
alpha = c(x2, a, b)

while (cond(x1, x2, eps) != 1):
    
    while abs(tGrad(x2, c(x2, a, b))) >= eps:
        
        if tGrad(x2, c(x2, a, b)) < 0:
            a = c(x2, a, b)
        if tGrad(x2, c(x2, a, b)) > 0:
            b = c(x2, a, b)        
    alpha = c(x2, a, b)

    for i in range(n):
        x1[i] = x2[i]
        x2[i] = x1[i] - alpha * MNeuton(x1)[i]
    k += 1
    
print(k, x2)
            

    
        
        
        
