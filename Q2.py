import numpy as np
import matplotlib.pyplot as plt

import endsem_library


def legendre(i ,X):
    if i == 0:
        return 1
    elif i == 1:
        return X
    elif i == 2:
        return (3 * X**2 - 1) / 2
    elif i ==3:
        return (5 * X**3 - 3*X) / 2
    elif i == 4:
        return (35 * X**4 - 30 * X**2 + 3) / 8


def func_leg(x, params):
    return sum([params[i] * legendre(i, x) for i in range(len(params))])


def polyfit_leg(x, y, sig, deg):
    n = len(x)
    A = [[0 for j in range(deg + 1)] for i in range(deg + 1)]
    b = [0 for i in range(deg + 1)]

    for i in range(n):
        for j in range(deg + 1):
            b[j] += legendre(j, x[i]) * y[i] / (sig[i] ** 2)
            for k in range(deg + 1):
                A[j][k] += (legendre(j, x[i]))*(legendre(k, x[i])) / (sig[i] ** 2)

    lower, upper = endsem_library.lu_decomposition(A)
    params = endsem_library.forward_backward_substitution(lower, upper, b)

    return params

mat = endsem_library.read_matrix("esem4fit.txt")
x, y = [], []
for i in range(len(mat)):
    x.append(mat[i][0])
    y.append(mat[i][1])

plt.scatter(x, y, label='datapoints')
plt.xlabel("x-axis")
plt.ylabel("y-axis")
plt.legend()
plt.show()
print("The scatter plot shows that the datapoints have three maxima+minima (in range -1 to 1). So, ideal choice of degree is 4")

params = polyfit_leg(x, y, sig=[1 for i in range(len(x))], deg=4)
print("Parameters:")
print(params)


y_n1 = []
for i in range(len(x)):
    y_n1.append(func_leg(x[i], params))
plt.scatter(x, y, s=1, label='datapoints')
plt.plot(x,y_n1, label='fitted curve with legendre basis')
plt.xlabel("x-axis")
plt.ylabel("y-axis")
plt.legend()
plt.show()


'''
The scatter plot shows that the datapoints have three maxima+minima (in range -1 to 1). So, ideal choice of degree is 4
Parameters:
[0.06965779687186321, 0.0036240203429268145, -0.012082580199521747, 0.011426217647052553, 0.11049235140900247]
'''

