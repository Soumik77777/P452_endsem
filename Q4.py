import numpy as np
import matplotlib.pyplot as plt

import endsem_library


def modi_func(input_func, y, a, b):
    y2 = ((a+b)/2) + ((b-a)*y/2)
    return input_func(y2) * ((b-a)/2)


def leggauss(deg):
    if deg == 4:
        x = [0.861136311, 0.339981043, -0.339981043, -0.861136311]
        w = [0.347854845, 0.652145154, 0.652145154, 0.347854845]
        return x, w

    elif deg == 5:
        x = [0.906179845, 0.538469310,  0.0        , -0.538469310,  -0.906179845]
        w = [0.236926885, 0.478628670, 0.568888889, 0.478628670, 0.236926885]
        return x, w
    elif deg == 6:
        x = [0.932469514, 0.661209386, 0.238619186,  -0.238619186,  -0.661209386, -0.932469514]
        w = [0.171324492, 0.360761573, 0.467913934, 0.467913934, 0.360761573, 0.171324492]
        return x, w


def integration_leg_gauss(input_func, deg, l, u):
    x, w = leggauss(deg)
    sum = 0
    for i in range(len(x)):
        sum += w[i] * modi_func(func, x[i], l, u)
    return sum


def func(x):
    return 1/(np.sqrt(1+(x**2)))

print("Potential at unit distance= ")
print("deg= 4, value= " +str(integration_leg_gauss(func, 4, -1, 1)))
print("deg= 5, value= " +str(integration_leg_gauss(func, 5, -1, 1)))
print("deg= 6, value= " +str(integration_leg_gauss(func, 6, -1, 1)))


'''
Potential at unit distance= 
deg= 4, value= 1.7620541789046658
deg= 5, value= 1.7628552954010728
deg= 6, value= 1.7627300484997592

# COMMENTS:
if we compare the three solutions with degree= 4, 5, 6; we can see that results with 5 and 6 degree is more accurate than degree 4.
'''
