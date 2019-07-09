import numpy as np
import random as rd
from scipy.spatial import distance


class Molecule:
    def __init__(self, coordinates):
        self.coordinates = coordinates


def select_process(constant_list):
    """
    CHOOSE A PATHWAY
    constant_list: rate constants list
    return: chosen path index
    """
    r = np.sum(constant_list) * rd.random()
    list_ = np.where(r > np.cumsum(constant_list))
    return len(list_[0])


def time_advance(constant_list):
    """
    TIME ADVANCE
    constant_list: list of rate constants
    return: advanced time
    """
    r = rd.random()
    return -np.log(r)/np.sum(constant_list)


def neighbourhood(center, radius, positions):
    """
    :param center: class
    :param radius: float
    :param positions: list
    :return: Elements in positions closer than a distance radius from center
    """
    center_position = np.array(center.coordinates)

    neighbours = []
    for element in positions:
        element_position = np.array(element.coordinates)

        if 0 < distance.euclidean(element_position, center_position) <= radius:
            neighbours.append(element)

    return neighbours


def particular_rates(center, directional_rates, neighbours):
    """
    :param center:  class
    :param directional_rates:   array
    :param neighbours: list of classes
    :return:  Hopping rate to each element of neighbours
    """
    center_position = np.array(center.coordinates)

    neighbour_rates = []

    for element in neighbours:
        element_position = np.array(element.coordinates)
        direction = (element_position - center_position) / distance.euclidean(element_position, center_position)

        element_rate = np.inner(direction, directional_rates)
        neighbour_rates.append(element_rate)

    return neighbour_rates



"""
Proposam una nova implementació del mètode de KMC en aquest cas tractant cada molècula
com una classe que es podrà extendre. Per començar només caracteritzada amb les seves coordenades (2D)
i sense cap mètode definit. Al no tenir més informació que les posicions prendrem els 'rates' com a
dades conegudes.
Disposarem les molècules en els nusos d'una xarxa cristal·lina (N, N).
"""

# Lattice generation
molecules = []
for i in range(1, 21):
    for j in range(1, 21):
        molecule = Molecule([i, j])
        molecules.append(molecule)

"""
Simulam la propagació d'una excitació. Aquesta pot decaure a l'estat fonamental o passar a una altra
molècula. Aquesta probabilitat de transició dependrà de la posició de la molècula. Vendrà donada
per l'expressió:
    k = (vect_k * r_director) / R
On r_director és el vector unitari en la direcció que uneix les dues molecules, vect_k és un vector
que modela les probabilitats de transició en cada direcció cartesiana (X, Y) privilegiant Y>0 i
R és la distància entre les molecules implicades.
Amb aquest disseny pensam que podrem passar al cas de distribució no ordenada més fàcilment.
"""

# Rate vector
rate_vector = np.array([1/7, 1/5])

decay_rate = 1/3

maximus_lifetime = 25

effective_radius = 1.05


"""
Prenem 10000 casos possibles de simulació. Ens quedarem amb les posicions de 
les molècules on ha arribat l'excitació i en farem l'estadística.
"""

exciton_paths = []
shift_list = []
time_list = []

for i in range(1000):
    path = []

    molecule = molecules[rd.randint(0, len(molecules) - 1)]
    path.append(molecule.coordinates)

    time = 0.0
    while time <= maximus_lifetime:

        neighbours = neighbourhood(molecule, effective_radius, molecules)
        neighbour_rates = (particular_rates(molecule, rate_vector, neighbours))
        neighbour_rates.insert(0, decay_rate)
        time = time + time_advance(neighbour_rates)

        process_index = select_process(neighbour_rates)
        if process_index == 0:
            break

        else:
            molecule = neighbours[process_index-1]
            path.append(molecule.coordinates)

    shift = distance.euclidean(np.array(path[0]), np.array(path[-1]))
    shift_list.append(shift)
    time_list.append(time)
