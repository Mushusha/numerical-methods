import numpy as np

def k22 (x, y):
    return x**2 + y**2 + 1

def func (v, i, j):
    return Ai(v, i, j) - b(i * h, j * h)
    
def b (y, x):
    #return -2*((y**2 - y)*(3*x**2 + y**2 - x + 1) + x**2 - x)
    return np.pi*np.sin(np.pi*y)*(-2*x*np.cos(np.pi*x)+(2+x**2+y**2)*np.pi*np.sin(np.pi*x))

def u (x, y):
    #return x * (x - 1) * y * (y-1)
    return np.sin(np.pi*x)*np.sin(np.pi*y)
k11 = 1
N = 10
l = 1
h = l / N
eps = 0.001
n = 1
def lambda_xx (v, i, j):
    return ((v[i + 1][j] - 2 * v[i][j] + v[i - 1][j])) / (h**2)

def lambda_yy (v, i, j):
    return (v[i][j + 1] * dRk22(i, j) - v[i][j] * (dRk22(i, j) + dLk22(i, j)) + v[i][j - 1] * dLk22(i, j)) / (h**2)

def dRk22 (i, j):
    return (k22(i * h, j * h) + k22(i * h, (j + 1) * h)) / 2

def dLk22 (i, j):
    return (k22(i * h, j * h) + k22(i * h, (j - 1) * h)) / 2

def Ai (v, i, j):
    #return - (diffLx(k11 * diffRx(v, i, j), i, j) + diffRx(k11 * diffLx(v, i, j), i, j) + diffLy(k22(i, j) * diffRy(v, i, j), i, j) + diffRy(k22(i, j) * diffLy(v, i, j), i, j)) / 2
    return -(lambda_xx(v, i, j) + lambda_yy(v, i, j))
    
def A (v):
    A_vn = np.zeros((N + 1, N + 1))
    for i in range(1, N):
        for j in range(1, N):
            A_vn[i][j] = -(lambda_xx(v, i, j) + lambda_yy(v, i, j))
            #A_vn[i][j] = - (diffLx(k11 * diffRx(v, i, j), i, j) + diffRx(k11 * diffLx(v, i, j), i, j) + diffLy(k22(i, j) * diffRy(v, i, j), i, j) + diffRy(k22(i, j) * diffLy(v, i, j), i, j)) / 2
    return A_vn


def alpha (j):
    c1 = 1
    c2 = 3
    return 1 / ((c1 + c2) / 2 + (c2 - c1) / 2 * np.cos(np.pi * (2 * j - 1) / (2 * n)))

def fill_muk2(v, muk2):
    for k2 in range(1, N):
        for i in range(1, N):
            muk2[i][k2] = mu_k2(v, i, k2)

def mu_k2(v, i, k2):
    sum = 0
    for j in range(1, N):
        sum += func(v, i, j) * np.sin(k2 * np.pi * j / N)
    return sum

def mu_k1k2(v, k1, k2, muk2):
    sum = 0
    for i in range(1, N):
        sum += muk2[i][k2] * np.sin(k1 * np.pi * i / N)
    return sum

def v_k2(v, i, k2, muk2):
    sum = 0
    for k1 in range(1, N):
        sum += mu_k1k2(v, k1, k2, muk2) / sobstv_ch(k1, k2) * np.sin(k1 * np.pi * i / N)
    return sum

def sobstv_ch(k1, k2):
    return (4 / h**2) * np.sin(k1 * np.pi * h /(2 * l))**2 + (4 / h**2) * np.sin(k2 * np.pi * h / (2 * l))**2

def v_(v, i, j, muk2):
    sum = 0
    for k2 in range(1, N):
        sum += v_k2(v, i, k2, muk2) * np.sin(k2 * np.pi * j / N)
    return (4 / (N**2)) * sum


v0 = np.zeros((N + 1, N + 1))
v = v0

muk2=np.zeros((N, N))

x = np.linspace(0, l, N + 1)
y = np.linspace(0, l, N + 1)
X, Y = np.meshgrid(x, y)

B = np.transpose(b(X, Y))
B[:, [0, -1]] = 0
B[[0, -1], :] = 0

U = np.transpose(u(X, Y))
err = np.max(np.abs(A(v) - B))
for k in range (n):
    if(err > eps):
        fill_muk2(v, muk2)
        y = np.zeros((N + 1,N + 1))
        for i in range(1, N):
            for j in range(1, N):
                y[i][j] = v_(v, i, j, muk2)       
        print(f'norma rn:{err} ')
        v -= alpha(k) * y
        err = np.max(np.abs(A(v) - B))
    else:
        break

while(err > eps):
    for k in range (n):
        fill_muk2(v, muk2)
        y = np.zeros((N + 1,N + 1))
        for i in range(1, N):
            for j in range(1, N):
                y[i][j] = v_(v, i, j, muk2)       
        print(f'norma rn:{err} ')
    v -= alpha(k) * y
    err = np.max(np.abs(A(v) - B))

#print(f'u на сетке:\n {U}')
#print(f'vn :\n {vn}')
print(f'||vn-u||={np.max(np.abs(v-u(X,Y)))}')
