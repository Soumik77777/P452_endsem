import numpy as np
import copy
import matplotlib.pyplot as plt

import endsem_library


def init_func(x):
    return 20 * abs(np.sin(x * np.pi))

def explicitPDE(inispace, a, b, k, dt, dx, tsteps=200):
    N = len(inispace)
    alpha = (k * dt) / (dx ** 2)
    for t in range(1, tsteps + 1):
        y_next = np.zeros(N)
        for i in range(1, N - 1):
            y_next[i] = inispace[i] + alpha * (inispace[i + 1] + inispace[i - 1] - 2 * inispace[i])
            y_next[0] = a*(t * dt)
            y_next[-1] = b*(t * dt)
        inispace = copy.deepcopy(y_next)

    return y_next

x = np.linspace(0, 2, 21)
init_T = [init_func(x[i]) for i in range(len(x))]
plt.plot(x, init_T, label='time step= 0')
t_list = [10, 20, 50, 100, 200, 500]
values = []
for i in range(len(t_list)):
    values.append(explicitPDE(init_T, 0, 0, 0.08, 0.008, 0.1, tsteps=t_list[i]))
    plt.plot(x, values[-1], label='time steps= '+str(t_list[i]))

plt.xlabel("length axis")
plt.ylabel("Temp profile (celcius)")
plt.legend()
plt.grid()
plt.show()


'''
The temperature profile at t=0 is a curve for modulus of sin curve. It is zero at both ends and centre.
The boundary points stays zero at all points. While all the centre points evolves as a function of
two different diffusion equations, one from 0 to l/2, another one from l/2 to l.
Both the lobes at t=0 could be evolved separately while the centre points are given as a composite of both of that.
'''
