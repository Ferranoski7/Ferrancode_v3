import matplotlib.pyplot as plt
from file_ed.read import read_txt
from electronic_density import ElectronicDensity
from molecule import Molecule
from temp import TemporalFunction
from symmetry_operations import SymOp
import numpy as np
from wfnsympy import WfnSympy

benzene = Molecule(['H','H'],[[-1,0,0],[1,0,0]],'Hydrogen')
sym_l = benzene.elements
coord_a = benzene.coordinates
name = benzene.name

Temp=SymOp(benzene)

soev = [SymOp(benzene).xcy(0, 6), SymOp(benzene).xcy(2, 6), SymOp(benzene).xcy(2, 3), SymOp(benzene).xcy(1, 2)]
selfsym = SymOp(benzene).analytic_selfsym()
benzene = ElectronicDensity(benzene)

basis = {
    'name': name,
    'primitive_type': 'gaussian',
}

alpha_mo_coeff = [[]]
atoms_data = []
for atom in range(len(sym_l)):
    shells_data = []
    element = sym_l[atom]
    p_exp = []
    con_coef = []
    p_con = []
    for shell in range(benzene.num_weights(element)):
        norm = 1 / np.sqrt(selfsym)

        p_exp.append(benzene.atomic_values[element][2 * (shell + 1)])
        con_coef.append((2 * benzene.atomic_values[element][2 * (shell + 1)] / np.pi) ** (3 / 2) * \
                        benzene.atomic_values[element][2 * (shell + 1) + 1]*benzene.atomic_values[element][0])
        p_con.append(0)

   # p_exp=[0.5]
    #con_coef=[1*(1 / np.pi) ** ( 3 / 2)]
    #p_con=[0]
    shells_data.append({
        'shell_type': 's',
        'p_exponents': p_exp,
        'con_coefficients': con_coef,
        'p_con_coefficients': p_con,
    })

    atoms_data.append({'shells': shells_data,
                       'symbol': element,
                       'atomic_number': benzene.atomic_values[element][0]})

    alpha_mo_coeff[0].append(norm)

basis['atoms'] = atoms_data

print(basis)
print(alpha_mo_coeff)
symmetry = WfnSympy(coord_a, sym_l, basis, alpha_mo_coeff, group='c2', axis=[0, 0, 1],center=[0,0,0])

symmetry.print_overlap_mo_alpha()

print(symmetry.self_similarity)

print(selfsym)


# Aquesta part del codi nomes era per comprovar que no donaven el mateix les autosemblances.
a = 2* 0.5

normi = 1 * benzene.atomic_values[element][0] * (a / np.pi) ** ( 3 / 2)

coef = (np.pi / (a + a)) ** (3 / 2)
expn = -1.0 * a * a / (a + a)
r2 = 0
ex = np.exp(expn * r2)

selfsim = normi * normi * coef * ex









exit()







































llista = ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']

name = 'At'
element = [name]
coord_a = [[0, 0, 0]]
element = Molecule(element, coord_a, name)
temp_element = TemporalFunction(element)

# wv_new = temp_element.gain_weights()

rho, x = temp_element.elden.linear_density(steps=20001)
# valence, _ = temp_element.linear_densityt(w=wv_new)

plt.figure('Separation_Core_and_Valence_' + name)
plt.plot(x, rho, label=' linear density')
# plt.plot(x, valence, label='New Valence density')
plt.legend()
plt.xlabel('r(A)')
plt.ylabel('Radial distribution')
plt.show()

for atom in range(0):
    # for atom in range(len(llista)):
    element = [llista[atom]]
    coord_a = [[0, 0, 0]]
    element = Molecule(element, coord_a, llista[atom])
    temp_element = TemporalFunction(element)

    wv_new = temp_element.gain_weights()
