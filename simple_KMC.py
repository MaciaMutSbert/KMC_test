#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 19:16:28 2019

@author: maciamutsbert

KINETIC MONTE CARLO SIMPLE EXAMPLE
"""

import numpy as np
import matplotlib.pyplot as plt
import random as rd
from matplotlib.ticker import NullFormatter

"""
DEFINING FUNCTIONS
"""


# Distribution function of firsts scape:
def probability(ki, t):
    return ki * np.exp(-ki * t)


def select_process(constant_list):
    """
    CHOOSE A PATHWAY
    constant_list: rate constants list
    return: chosen path index
    """
    r = np.sum(constant_list) * rd.random()
    list_ = np.where(r > np.cumsum(constant_list))
    return len(list_[0])+1


def time_advance(constant_list):
    """
    TIME ADVANCE
    constant_list: list of rate constants
    return: advanced time
    """
    r = rd.random()
    return -np.log(r)/np.sum(constant_list)


def count_and_remove_0(any_list):
    number_of_0 = 0
    index_list = []
    for i, li in enumerate(any_list):
        if li == 0:
            number_of_0 = number_of_0+1
            index_list.append(i)

    non_zero_list = np.delete(any_list, index_list)
    return number_of_0, non_zero_list
            

"""
-----------KMC----------
Suposam que tenim una xarxa (2D) de molecules, totes idèntiques entre sí.
Per major simplicitat la suposarem infinita, i prendrem uns eixos de coordenades.
El transport serà diferent segons l'eix escollit.
Se'n excita una qualsevol (la suposarem a (0,0)) al seu primer estat excitat. 
Hi ha 5 camins possibles:
    a) L'electró cau a l'estat fonamental.
    b) L'excitació es mou segons l'eix X (cap a x>0 o cap a x<0)
    c) L'excitació es mou segons l'eix Y (cap a y>0 o cap a y<0)
"""
# Rates
rate_1 = 1/6  # ground stage decay rate
rate_2 = 1/5  # x>0 transport rate
rate_3 = 1/5  # x<0 transport rate
rate_4 = 1/7  # y>0 transport rate
rate_5 = 1/7  # y<0 transport rate

rate_constants = np.array([rate_1, rate_2, rate_3, rate_4, rate_5])
t_max = 50.0


"""
Prenem 10000 casos possibles de simulació. Ens quedarem amb les posicions de 
les molècules on ha arribat l'excitació i en farem l'estadística.
"""

x_list = []
y_list = []
time_list = []
for i in range(10000):
    
    x_shift = 0
    y_shift = 0
    time = 0.0
    while time < t_max:
        process_index = select_process(rate_constants)
        if process_index == 1:
            break
        
        elif process_index == 2:
            x_shift = x_shift+1
            
        elif process_index == 3:
            x_shift = x_shift-1
            
        elif process_index == 4:
            y_shift = y_shift+1
            
        elif process_index == 5:
            y_shift = y_shift-1
            
        time = time + time_advance(rate_constants)
        
    x_list.append(x_shift)
    y_list.append(y_shift)
    time_list.append(time)

x_list = np.array(x_list)
y_list = np.array(y_list)
time_list = np.array(time_list)


"""
                Stadistic study
"""
# 1. Maximum displacements
x_absolute_shift = np.abs(x_list)
y_absolute_shift = np.abs(y_list)

print('Maximum displacement along X axis:', x_absolute_shift.max())
print('Maximum displacement along Y axis:', y_absolute_shift.max())


# 2. Average displacement

print('Average x displacement: ', np.average(x_list))
print('Average y displacement: ', np.average(y_list))

# 3. Standard deviation
x_deviation = np.sqrt(np.average(x_list**2) - np.average(x_list)**2)
y_deviation = np.sqrt(np.average(y_list**2) - np.average(y_list)**2)

print('X standard deviation: ', x_deviation)
print('Y standard deviation: ', y_deviation)

# 4. Non-propagate excitons
radius = np.sqrt(x_list**2 + y_list**2)

print('Non-propagate excitons: ', count_and_remove_0(radius)[0])

# 5. Average lifetime of the exciton and its standard deviation

print('Average lifetime of the exciton: ', np.average(time_list))

time_deviation = np.sqrt(np.average(time_list**2)-np.average(time_list)**2)

"""
                        Coeficient de difusió
Podem definir un coeficient de difusió a partir de la distància mitjana 
recorreguda i el temps de vida mitjà de l'excitó en virtud de la relació:
    L=sqrt(D*t)
    on D és l'esmentat coeficient de difusió
        D=(L**2)/t
"""

diffusion_constant = (np.average(radius)**2)/np.average(time_list)

print('Difussion coeficient of the lattice: ', diffusion_constant)


"""
                2D histogram
"""
# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left + width + 0.02

hist_2d_rect = [left, bottom, width, height]
x_hist_rect = [left, bottom_h, width, 0.2]
y_hist_rect = [left_h, bottom, 0.2, height]

# start with a rectangular Figure
plt.figure(1, figsize=(8, 8))

ax_hist2d = plt.axes(hist_2d_rect)
ax_hist_x = plt.axes(x_hist_rect)
ax_hist_y = plt.axes(y_hist_rect)

# no labels
ax_hist_x.xaxis.set_major_formatter(NullFormatter())
ax_hist_y.yaxis.set_major_formatter(NullFormatter())


# now determine nice limits by hand:
binwidth = 0.45
xymax = np.max([np.max(np.fabs(x_list)), np.max(np.fabs(y_list))])
lim = (int(xymax/binwidth) + 1) * binwidth

ax_hist2d.set_xlim((-lim, lim))
ax_hist2d.set_ylim((-lim, lim))

bins = np.arange(-lim, lim + binwidth, binwidth)
ax_hist_x.hist(x_list, bins=bins, color='red')
ax_hist_x.set_xlim(ax_hist2d.get_xlim())
ax_hist_x.set_title('Exciton final positions histogram')
ax_hist_x.set_ylabel('X frequency')

ax_hist_y.hist(y_list, bins=bins, orientation='horizontal', color='red')
ax_hist_y.set_ylim(ax_hist2d.get_ylim())
ax_hist_y.set_xlabel('Y frequency')

# plot:
ax_hist2d.hist2d(x_list, y_list, bins=bins, cmap='Reds')
ax_hist2d.set_xlabel('X')
ax_hist2d.set_ylabel('Y')
plt.show()
plt.tight_layout()


"""
                    Final positions plot
"""
plt.plot(x_list, y_list, 'ro', label='Position')
plt.xlabel('X')
plt.ylabel('Y')
plt.legend()
plt.title('Excitons final positions')
plt.show()