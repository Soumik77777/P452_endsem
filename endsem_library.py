import math
import csv
import numpy as np
import scipy
import matplotlib.pyplot as plt



### Matrix Operations

## Read Matrix
def read_matrix (txt):
    with open(txt, 'r') as a:
        matrix = [[float(num) for num in row.split(' ')] for row in a]
    return matrix


## Matrix multiply
def mat_multiply(mat_a, mat_b):
    if len(mat_a[0])!=len(mat_b):
        print("Recheck the matrix shapes.")
        return 0
    else:
        mat_c = np.zeros((len(mat_a), len(mat_b[0])))
        for i in range(len(mat_a)):
            for j in range(len(mat_b[0])):
                sum = 0
                for k in range(len(mat_a[0])):
                    sum += mat_a[i][k] * mat_b[k][j]
                mat_c[i,j] = sum
        return mat_c




## Partial Pivot
def partial_pivot (m, v):
    n = len(m)
    for i in range (n-1):
        if m[i][i] ==0:
            for j in range (i+1,n):
                if abs(m[j][i]) > abs(m[i][i]):
                    m[i], m[j] = m[j], m[i]
                    v[i], v[j] = v[j], v[i]
    return m,v



## LU decomposition
def lu_decomposition(matrix):
    n = len(matrix)
    upper_mat = [[0 for i in range(n)] for j in range(n)]
    lower_mat = [[0 for i in range(n)] for j in range(n)]

    for i in range(n):
        for j in range(i, n): #calculating upper matrix
            sum = 0
            for k in range(i):
                sum += (lower_mat[i][k] * upper_mat[k][j])
            upper_mat[i][j] = matrix[i][j] - sum

        for j in range(i, n): #calculating lower matrix
            if (i == j):
                lower_mat[i][i] = 1     ## Doolittle
            else:
                sum = 0
                for k in range(i):
                    sum += (lower_mat[j][k] * upper_mat[k][i])

                lower_mat[j][i] = ((matrix[j][i] - sum) / upper_mat[i][i])

    return lower_mat, upper_mat



## Cholesky decomposition
def cholesky_decompose(matrix):
    n = len(matrix)
    lower = np.zeros_like(matrix)

    lower[0, 0] = np.sqrt(matrix[0][0])

    for j in range(1, n):
        lower[j, 0] = matrix[j][0] / lower[0, 0]

    for i in range(1, n):
        sum = 0.0
        for p in range(i):
            sum = sum + lower[i, p] ** 2
        lower[i, i] = np.sqrt(matrix[i][i] - sum)

        for j in range(i + 1, n):
            sum = 0.0
            for p in range(i):
                sum = sum + lower[i, p] * lower[j, p]
            lower[j, i] = matrix[j][i] - sum / lower[i, i]

    return lower, np.matrix.transpose(lower)



## forward-backward substitution
def forward_backward_substitution (lower_mat, upper_mat, vector):
    '''
    If we have LUx=B,
    first we solve Ly=B, then Ux=y
    '''
    n = len(lower_mat)
    # forward-substitution
    y = [0] * n
    for i in range(n):
        sum = 0
        for j in range(i):
            sum += lower_mat[i][j] * y[j]

        y[i] = vector[i] - sum

    #backward-substitution
    x = [0] * n
    for i in reversed(range(n)):
        sum = 0
        for j in range(i + 1, n):
            sum+= upper_mat[i][j] * x[j]
        x[i] = (y[i] - sum)/ upper_mat[i][i]

    return x



def inverse_lu(matrix):
    n = len(matrix)
    lower, upper = lu_decomposition(matrix)
    inv_mat = np.zeros((n,n))

    for i in range(n):
        vector = [0 for i2 in range(n)]
        vector[i] = 1
        x = forward_backward_substitution(lower, upper, vector)
        for j in range(n):
            inv_mat[j, i] = x[j]
    return inv_mat





# ---------------------------------------------------------------------------------------


### Curve Fitting

## chi-sqr
def chi_square(x_list, y_list, exp_data):
    chi_sqr = 0
    for i in range(len(x_list)):
        chi_sqr += (y_list[i] - exp_data[i]) ** 2 / exp_data[i]
    return chi_sqr



def polyfit(x, y, sig, deg):
    n = len(x)
    A = [[0 for j in range(deg + 1)] for i in range(deg + 1)]
    b = [0 for i in range(deg + 1)]

    for i in range(n):
        for j in range(deg + 1):
            b[j] += (x[i] ** j) * y[i] / (sig[i] ** 2)
            for k in range(deg + 1):
                A[j][k] += (x[i] ** (j + k)) / (sig[i] ** 2)

    lower, upper = lu_decomposition(A)
    params = forward_backward_substitution(lower, upper, b)

    return params

'''
mat = read_matrix("esem4fit.txt")
x, y = [], []
for i in range(len(mat)):
    x.append(mat[i][0])
    y.append(mat[i][1])
params1 = polyfit(x, y, sig=[1 for i in range(len(x))], deg=4)
'''



def chebyshev(i ,X):
    if i == 0:
        return 1
    elif i == 1:
        return 2 * X - 1
    elif i == 2:
        return 8 * X ** 2 - 8 * X + 1
    elif i ==3:
        return 32 * X ** 3 - 48 * X ** 2 + 18 * X - 1
    elif i == 4:
        return 128*X**4 - 256*X**3 + 160*X**2 - 32*X + 1


def polyfit_chev(x, y, sig, deg):
    n = len(x)
    A = [[0 for j in range(deg + 1)] for i in range(deg + 1)]
    b = [0 for i in range(deg + 1)]

    for i in range(n):
        for j in range(deg + 1):
            b[j] += chebyshev(j, x[i]) * y[i] / (sig[i] ** 2)
            for k in range(deg + 1):
                A[j][k] += (chebyshev(j, x[i]))*(chebyshev(k, x[i])) / (sig[i] ** 2)

    lower, upper = lu_decomposition(A)
    params = forward_backward_substitution(lower, upper, b)

    return params

'''
params2 = polyfit_chev(x, y, sig=[1 for i in range(len(x))], deg=4)
def func(x, params):
    sum = 0
    for i in range(len(params)):
        sum += params[i] * x **i
    return sum
def func2(x, params):
    return sum([params[i] * chebyshev(i, x) for i in range(len(params))])


y_n1, y_n2 = [], []
for i in range(len(x)):
    y_n1.append(func(x[i], params1))
    y_n2.append(func2(x[i], params2))
plt.scatter(x, y, s=1)
plt.plot(x,y_n1)
plt.plot(x, y_n2)
plt.show()
'''




### Random number

## pseudo rng
def pRNG(seed, a, m, n):
    x = [seed]
    for i in range(n):
        x.append((a*x[i])%m)
    x = [i/m for i in x]
    return x


## random walk
def random_walk_2d(N):                                  ##Function for random walk
    x, y = 0, 0
    x_positions, y_positions = [x], [y]
    for i in range(N):
        z = random.uniform(0, 2 * math.pi)              ##generating random angle in range (0,2*pi)
        x += math.cos(z)                                ##displacement in x is cos component
        y += math.sin(z)                                ##displacement in y is sin component
        x_positions.append(x)
        y_positions.append(y)
    position = [x_positions,y_positions]
    return position, x, y                             ##returning a list of position after each step and last coordinate






### Gaussian quadrature

## func

## modified func

## cheb1

## cheb2

## cheb3
# W = np.sqrt(1+x/1-x)



## cheb4

## legendre
# W = 1
# xi = cos(pi * (0.5+i)/N )
# wi = 2 / ( () * ())
# np.polynomials.legendre.leggauss
## integration





### pde







### coupled ODE

## RK4

## coupled RK4















