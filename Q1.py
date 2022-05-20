import numpy as np
import matplotlib.pyplot as plt

import endsem_library


## pseudo random number generator with linear congruent method, normalized from 0 to 1
def pRNG(seed, a, m, N):
    x = [seed]
    for i in range(N):
        x.append((a*x[i])%m)
    x = [i/m for i in x]
    return x


## code for random walk
def random_walk_2d(seed, a, m, N):                                  ##Function for random walk
    x, y = 0, 0
    x_positions, y_positions = [x], [y]
    random_list = pRNG(seed, a, m, N-1)
    random_angles = [random_list[i] * 2*np.pi for i in range(len(random_list))]
    for i in range(N):
        x += np.cos(random_angles[i])                                ##displacement in x is cos component
        y += np.sin(random_angles[i])                                ##displacement in y is sin component
        x_positions.append(x)
        y_positions.append(y)

    return x_positions, y_positions                         ##returning a list of x,y coordinates after each step



## generate and plotting random walk of 200 steps
seeds = [119, 4096, 8312, 11939, 13789]
for i in range(len(seeds)):
    random1 = random_walk_2d(seeds[i], 572, 16381, 200)     # seed chosen manually
    plt.plot(random1[0], random1[1], label='Random walk number '+str(i+1))
    plt.xlabel("X-direction")
    plt.ylabel("Y-direction")

plt.title('Random walk in 2D for 200 steps')
plt.legend(shadow=True)
plt.grid()
plt.show()


## relation btwn RMS and N
R_avg = []
step_list = [200, 500, 1000, 1500, 2000]
for j in range(len(step_list)):
    list_R = []
    list_seed = np.linspace(1, 16381, 100)
    for i in range(100):
        x_list, y_list = random_walk_2d(list_seed[i], 572, 16381, step_list[j])
        R = np.sqrt(x_list[-1] ** 2 + y_list[-1] ** 2)  ##for rms value
        list_R.append(R)
    R_avg.append(np.mean(list_R))

root_N = [np.sqrt(step_list[i]) for i in range(len(R_avg))]
plt.plot(root_N, R_avg)
plt.xlabel("Square root of N")
plt.ylabel("rms of Radial Distance for each N")
plt.title('Plot of R_rms vs sqrt of N')
plt.grid()
plt.show()

print("Average RMS for 100 different random walks= " + str(R_avg[0]))
print("Square root of 200= "+ str(np.sqrt(200)))





'''
Average RMS for 100 different random walks= 14.510593374439969
Square root of 200= 14.142135623730951

## additional comments
We can see that average of RMS is close to the root of N, 
Plotting a curve with changing N shows that the relation between root N and average RMS 
is a polynomial of order 1.
'''
