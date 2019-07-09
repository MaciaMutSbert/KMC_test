#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  6 20:20:19 2019

@author: maciamutsbert

POLYMER CHAIN MORPHOLOGY GENERATION BASED ON A SELF-AVOIDING RANDOM WALK
"""

import numpy as np
import random as rd


def available_positions(matrix, x, y):
    """
    Funció que donada la posicio (x, y) en una matriu M ens diu quantes i
    quines posicions valen 0 (monomers lliures). Per tant reb com arguments 
    una matriu i 2 índexs i retorna una llista amb les posicions lliures
    """
    possible_positions = []
    
    for [i, j] in [[x+1, y], [x-1, y], [x, y+1], [x, y-1]]:
        
        if matrix[i, j] == 0:
            possible_positions.append([i, j])

    return possible_positions


def initial_position(matrix, n):
    z = 0
    while z <= n:
        i, j = np.random.randint(0, 21, size=2)
        z += 1
        
        if matrix[i, j] == 0:
            break
    return [i, j], z


"""
Prenem una matriu (N, N) de 0s que serà la nostra xarxa de monòmers disponibles.
Per solventar el problema de les vores crearem una matriu (N+2, N+2) amb
valors -1 (que l'algorisme entendrà com a monòmers no disponibles).
"""

N = 20

monomer_lattice = np.zeros((N+2, N+2))
total_monomers = N**2

for i in (0, N+1):
    monomer_lattice[:, i] = -1
    monomer_lattice[i, :] = -1
     
    
#############################################

polymer_file = open('Polymers.txt', 'w')

polymers = []
counter = 0
row_num = 1


while counter <= total_monomers:
    
    # Escollim una posició inicial entre els monòmers lliures
    i, j = initial_position(monomer_lattice, total_monomers)[0]

    if initial_position(monomer_lattice, total_monomers)[1] == total_monomers:
        print('Lattice completed with polymers 2')
        break
        
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
    
        available_monomers = available_positions(monomer_lattice, i, j)
    
        if len(available_monomers) == 0:
            break
        
        else:
            # -------------------Nova posició---------------------
            [i, j] = available_monomers[rd.randint(0, len(available_monomers)-1)]
            
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
        break


polymer_file.close()



    
   
    
    



            