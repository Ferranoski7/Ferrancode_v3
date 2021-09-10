import matplotlib.pyplot as plt
from file_ed.read import read_txt
from electronic_density import ElectronicDensity
from molecule import Molecule
from temp import TemporalFunction
from symmetry_operations import SymOp
import numpy as np


llista = ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']


element = ['N']
coord_a = [[0, 0, 0]]
element=Molecule(element,coord_a,'N')
temp_element=TemporalFunction(element)

wv_new=temp_element.gain_weights()


rho, x = temp_element.elden.linear_valence_density()
valence, _ = temp_element.linear_densityt(w=wv_new)





plt.figure('Separation_Core_and_Valence_N')
plt.plot(x, rho, label='Valence linear density')
plt.plot(x, valence, label='New Valence density')
plt.legend()
plt.xlabel('r(A)')
plt.ylabel('Radial distribution')
plt.show()


for atom in range(0):
#for atom in range(len(llista)):
    element = [llista[atom]]
    coord_a = [[0, 0, 0]]
    element=Molecule(element,coord_a,llista[atom])
    temp_element=TemporalFunction(element)

    wv_new=temp_element.gain_weights()




