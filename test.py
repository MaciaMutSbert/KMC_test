import numpy as np
import matplotlib.pyplot as plt
import random as rd


# Distribution function of firts scape:
def probability(ki, t):
    return ki * np.exp(-ki * t)


def select_process(rate_constants):
    """
    CHOOSE A PATHWAY
    rate_constants: list of rate constants
    return: chosen path index
    """
    r = np.sum(rate_constants) * rd.random()
    index_list = np.where(r > np.cumsum(rate_constants))

    return len(index_list[0])+1


def time_advance(rate_constants):
    """
    TIME ADVANCE
    rate_constants: list of rate constants
    return: advanced time
    """
    r = rd.random()
    return -np.log(r)/np.sum(rate_constants)


"""
----------KMC----------
"""


# definim paramatres inicials
t_max = 100.0
rate_constants = np.array([1/5, 1/7, 1/4, 1/9])

time = 0.0
index_path = []
while time < t_max:
    process_index = select_process(rate_constants)
    index_path.append(process_index)
    time = time + time_advance(rate_constants)


# Probability distributions for each rate
time_range = np.arange(0, 50, 0.1)
for i, ki in enumerate(rate_constants):
    plt.plot(time_range, [probability(ki, t) for t in time_range], label='rate constant: {}'.format(i+1))

plt.legend()
plt.xlabel('Time')
plt.ylabel('Probabiliy')
plt.title('Probability distribution for each rate')
plt.show()


# Histogram
plt.hist(index_path)
plt.title('Histogram')
plt.xlabel('Chosen path')
plt.ylabel('Frequency')
plt.show()

