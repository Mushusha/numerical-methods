import numpy as np
import math

def func(x, y):
    #return np.array([y[1], y[0]])
    return np.array([y])

def eps(h):
    return e / (b - a)

a = 0
b = 1
x = np.array([a, a, a])
x = x.astype("float")
#y0 = np.array([0, 2])
y0 = np.array([1])
y = np.array([y0, y0, y0])
y = y.astype("float")
h = (b - a)
e = 10 ** (-3)
p = 3

def y_n(x, y, h):
    a = []
    a.append(y + (k1(x, y, h) + 3 * k3(x, y, h)[0]) / 4)
    a.append(y + (k1(x, y, h) + 4 * k3(x, y, h)[1] + k4(x, y, h)[0]) / 6)
    a.append(y + (k1(x, y, h) + 4 * k4(x, y, h)[1] + k5(x, y, h)) / 6)
    return a

def k1(x, y, h):
    return h * func(x, y) # 1 2 3

def k2(x, y, h):
    a = []
    a.append(h * func(x + h/3, y + k1(x, y, h) / 3)) # 1 3
    a.append(h * func(x + h/4, y + k1(x, y, h) / 4)) # 2
    return a

def k3(x, y, h):
    a = []
    a.append(h * func(x + 2*h/3, y + 2 * k2(x, y, h)[0] / 3)) # 1
    a.append(h * func(x + h/2, y + k2(x, y, h)[1] / 2)) # 2
    a.append(h * func(x + h/2, y + (k1(x, y, h) + k2(x, y, h)[0]) / 6)) # 3
    return a

def k4(x, y, h):
    a = []
    a.append(h * func(x + h, y + k1(x, y, h) - 2 * k2(x, y, h)[1] + 2 * k3(x, y, h)[1])) # 2
    a.append(h * func(x + h/2, y + (k1(x, y, h) + 3 * k3(x, y, h)[2]) / 8)) # 3
    return a

def k5(x, y, h):
    return h * func(x + h, y + (k1(x, y, h) - 3 * k3(x, y, h)[2]) / 2 + 2 * k4(x, y, h)[1]) # 3

def E(x, y, h):
    return (2 * k1(x, y, h) - 9 * k3(x, y, h)[2] + 8 * k4(x, y, h)[1] - k5(x, y, h)) / 30

def ro(x, y, h):
    a = []
    a.append(abs((y_n(x, y, h)[0] - y_n(x, y, h/2)[0]) / (1 - 2 ** (-p))))
    a.append(abs(y_n(x, y, h)[0] - y_n(x, y, h)[1]))
    a.append(abs(E(x, y, h)))
    return a

def check(a):
    for i in a:
        if i <= eps(h):
            return False
    return True

t = 1
for j in range (3):
    i = 0
    while b - x[j] > e:
        if check(ro(x[j], y[j], h)[j]) and t:
            h /= 2
        elif not(check(ro(x[j], y[j], h)[j])) and t:
            h *= 2
            t = 0
        y[j] = y_n(x[j], y[j], h)[j]
        x[j] += h
        i += 1
    print("Метод: ", j + 1, "\nЗначение в конце отрезка: ", y[j], "\nКоличество возвращений к правой части: ", i, sep = "", end = "\n\n")
