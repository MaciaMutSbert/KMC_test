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

"""
DEFINING FUNCTIONS
"""


# Distribution function of firts scape:
def probability(ki, t):
    return ki * np.exp(-ki * t)


def select_process(constants_list):
    """
    CHOOSE A PATHWAY
    constants_list: rate constants list
    return: chosen path index
    """
    r = np.sum(constants_list) * rd.random()
    list_ = np.where(r > np.cumsum(constants_list))
    return len(list_[0])+1


def time_advance(constants_list):
    """
    TIME ADVANCE
    constants_list: list of rate constants
    return: advanced time
    """
    r = rd.random()
    return -np.log(r)/np.sum(constants_list)


def count_0(l):
    
    number_of_0 = 0
    for i in l:
        if i == 0:
            number_of_0 = number_of_0+1
            
    return number_of_0
            

"""
                                KMC
Suposam que tenim una xarxa (2D) de molecules, totes idèntiques entre sí.
Per major simplicitat la suposarem infinita, i prendrem uns eixos de coordenades.
El trannsport serà diferent segons l'eix escollit.
Se'n excita una qualsevol (la suposarem a (0,0) al seu primer estat excitat. 
Hi ha 5 camins possibles:
    a) L'electró cau a l'estat fonamental.
    b) L'excitació es mou segons l'eix X (cap a x>0 o cap a x<0)
    c) L'excitació es mou segons l'eix Y (cap a y>0 o cap a y<0)
"""
# Rates
rate_1 = 1/3    # ground stage decay rate
rate_2 = 1/5    # x>0 transport rate
rate_3 = 1/5    # x<0 transport rate
rate_4 = 1/7    # y>0 transport rate
rate_5 = 1/7    # y<0 transport rate

rate_constants = np.array([rate_1,rate_2,rate_3,rate_4,rate_5])
t_max = 50.0


"""
Prenem 1000 casos possibles de simulació. Ens quedarem amb les posicions de 
les molècules on ha arribat l'excitació i en farem l'estadística.
"""

x_list = []
y_list = []
time_list = []
for i in range(1000):
    
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
                    Results
"""    
plt.plot(x_list, y_list, 'ro', label='Position')
plt.legend()
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Final positions')
plt.show()


"""
                Stadistic study
"""
# 1. Maximum displacements
x_absolute_shift = np.abs(x_list)
y_absolute_shift = np.abs(y_list)

print('Maximum displacement along X axis:', x_absolute_shift.max())
print('Maximum displacement along Y axis:', y_absolute_shift.max())


# 2. Average displacement
radius = np.sqrt(x_list**2 + y_list**2)

print('Average radius',np.average(radius))


# 3. Standard desviation
radius_deviation = np.sqrt(np.average(radius**2) - np.average(radius)**2)

print('Radius standard deviation:', radius_deviation)

# 4. Non-propagate excitons

print('Non-propagate excitons: ', count_0(radius))

# 5. Average lifetime of the exciton

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

difussion_constant = (np.average(radius)**2)/np.average(time_list)

print('Difussion coeficient of the lattice: ', difussion_constant)


