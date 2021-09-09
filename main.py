import matplotlib.pyplot as pltCancel changes
from file_ed.read import read_txt
from electronic_density import ElectronicDensity
from molecule import Molecule
from temp import TemporalFunction
from symmetry_operations import SymOp
import numpy as np


llista = ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']

for atom in range(len(llista)):
    element = [llista[atom]]
    coord_a = [[0, 0, 0]]

    element=Molecule(element,coord_a,llista[atom])
    temp_element=TemporalFunction(element)

    wv_new=temp_element.gain_weights()




