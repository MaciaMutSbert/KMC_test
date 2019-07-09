#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  7 20:19:40 2019

@author: maciamutsbert

POLYMER CHAIN MORPHOLOGY GENERATION BASED ON A SELF-AVOIDING RANDOM WALK
"""

import numpy as np
import random as rd


def available_positions(M, x, y, rate_constants):
    """
    Funció que donada la posicio (x, y) en una matriu M ens diu quantes i
    quines posicions valen 0 (monomers lliures) i la probabilitat de cada posició.
    Per tant reb com arguments una matriu, 2 índexs i una llista amb les 
    probabilitats de cada direcció.
    Retorna 2 llistes, una amb les posicions disponibles i una altre amb les 
    probabilitats d'aquestes. Estan ordenats igual.
    """
    positions = [[x+1, y], [x-1, y], [x, y+1], [x, y-1]]
    
    possible_positions = []
    available_rates = []

    for k in range(0, len(positions)):
        i = positions[k][0]
        j = positions[k][1]
        
        if M[i, j] == 0:
            possible_positions.append([i, j])
            available_rates.append(rate_constants[k])
            
    return possible_positions, available_rates


def select_process(rate_constants):
    """
    CHOOSE A PATHWAY
    rate_constants: rate constants list
    return: chosen path index
    """
    r = np.sum(rate_constants) * rd.random()
    list_ = np.where(r > np.cumsum(rate_constants))
    return len(list_[0])


def initial_position(matrix, n):
    z = 0
    while z <= n:
        i, j = np.random.randint(0, 21, size = 2)
        z += 1
        
        if matrix[i, j] == 0:
            break
    return [i, j]


"""
Programa que genera polimers a partir d'una xarxa (N, N) monòmers per mitjà
d'un algorisme "self-avoiding random walk" amb una direcció privilegiada.
Per la tria d'un dels camins possibles implementam l'algorisme KMC.
Per solventar el problema de les vores crearem una matriu (N+2, N+2) amb
valors -1 (que l'algorisme entendrà com a monòmers no disponibles).
"""


# Lattice generation
N = 20

monomer_lattice = np.zeros((N+2, N+2))
total_monomers = N**2

for i in (0, N+1):
    monomer_lattice[:, i] = -1
    monomer_lattice[i, :] = -1


"""
Definim les probabilitats d'incloure un monòmer segons la direcció, 
privilegiam la direcció y>0. Aquesta preferència es quantitza amb un paràmetre
mu en [0, 1]. Prendrem mu = 0,6. Els valors de les noves probabilitats seran:
    1+mu per la direcció privilegiada, y>0
    1-mu per les direccions paral·leles, x
    (1-mu/1+mu)**2 per la direcció contraria a la privilegiada, y<0
    
Per conveni en contruir la llista amb aquestes probabilitats seguirem el següent
ordre: [x+1_rate, x-1_rate, y+1_rate, y-1_rate]
"""

mu = 0.6
constant_rates = [1-mu, 1-mu, 1+mu, ((1-mu)/(1+mu))**2]

     
#############################################

polymer_file = open('Polymers_2.txt', 'w')

polymers = []
counter = 0
row_num = 1


while counter <= total_monomers:
    
    # Escollim una posició inicial entre els monòmers lliures
    i, j = initial_position(monomer_lattice, total_monomers)
        
    #######
    # Generam una cadena de monomers, un polimer
    #######
    
    chain = [[i, j]]
    monomer_lattice[i, j] = 1       # canviam de valor perquè deixi de comptar com a lliure
    counter += 1

    # Longitud, donada per una gaussiana de mitjana 6 i desviació 1.5
    chain_length = round(rd.gauss(6.0, 1.5))
    
    chained_monomers = 1
    while chained_monomers <= chain_length:
    
        available_monomers, available_rates = available_positions(monomer_lattice, i, j, constant_rates)
        
        if len(available_monomers) == 0:
            break
        
        else:
            """
            La nova posició s'ecollirà de manera aleatòria entre 
            les possibles amb l'algorisme de KMC
            """
            i, j = available_monomers[select_process(available_rates)]
            
            chain.append([i, j])
            monomer_lattice[i, j] = 1
            chained_monomers += 1
            counter += 1
    
    polymers.append(chain)
    
    # Guardam els resultats en un fitxer
    polymer_file.write('Polymer %2d \t' % row_num + str(chain) + '\n')
    row_num += 1
    
    if counter == total_monomers:
        print('Lattice completed with polymers 1')
        there_are_free_monomers = False


polymer_file.close()



    
   
    
    



            

