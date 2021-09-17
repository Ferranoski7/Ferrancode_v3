import matplotlib.pyplot as plt
from file_ed.read import read_txt
from file_ed.read import read_xyz
from electronic_density import ElectronicDensity
from molecule import Molecule
from temp import TemporalFunction
from symmetry_operations import SymOp
import numpy as np
from wfnsympy import WfnSympy
from rotation.rotate import rotation

benzene = read_xyz('coord')
#benzene = Molecule(['C', 'H'], [[-1, 0, 0], [1, 0, 0]], 'CH')
temp = SymOp(benzene)
print(temp.analytic_selfsym())

print('---------------------------------\n Our measures:\n',
      temp.xcy(6, 6, axis=[-1.0685,     -0.0537,      0.1921], cm=[0, 0, 0]),
      temp.xcy(1, 6, axis=[-1.0685,     -0.0537,      0.1921], cm=[0, 0, 0]),
      temp.xcy(1, 3, axis=[-1.0685,     -0.0537,      0.1921], cm=[0, 0, 0]),
      temp.xcy(1, 2, axis=[-1.0685,     -0.0537,      0.1921], cm=[0, 0, 0]))
print('---------------------------------\n Wyfnsym measures:\n')
sym_l = benzene.elements
coord_a = benzene.coordinates
name = benzene.name

Temp = SymOp(benzene)
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
        Element = Molecule([element], [[0, 0, 0]])
        selfsymi = SymOp(Element).analytic_cosymlib_integral()
        #print(selfsymi)
        expi = 2 * benzene.atomic_values[element][2 * (shell + 1)]
        coefi = benzene.atomic_values[element][2 * (shell + 1) + 1]
        nrmi = ((2*expi / np.pi) ** (3 / 4))

        p_exp.append(expi)
        con_coef.append(coefi * benzene.atomic_values[element][0] * nrmi)
        p_con.append(0.0)

    # p_exp=[0.5]
    # con_coef=[1*(1 / np.pi) ** ( 3 / 2)]
    # p_con=[0]
    shells_data.append({
        'shell_type': 's',
        'p_exponents': p_exp,
        'con_coefficients': con_coef,
        'p_con_coefficients': p_con,
    })

    atoms_data.append({'shells': shells_data,
                       'symbol': element,
                       'atomic_number': benzene.atomic_values[element][0]})
    alpha_mo_coeff[0].append(1.0* np.sqrt(selfsymi))

print(alpha_mo_coeff)
basis['atoms'] = atoms_data

# print(basis)

symmetry = WfnSympy(coord_a, sym_l, basis, alpha_mo_coeff, group='c6', axis=[-1.0685,     -0.0537,      0.1921],
                    center=[0, 0, 0])




print(symmetry.mo_SOEVs_a / symmetry.mo_SOEVs_a[0][0])

symmetry.print_overlap_mo_alpha()

exit()
from scipy import integrate



def g1(x, y, z, px, py, pz, exp):
    x, y, z, px, py, pz = x / 0.529177249, y / 0.529177249, z / 0.529177249, px / 0.529177249, py / 0.529177249, pz / 0.529177249

    return np.exp(-exp * ((x - px) ** 2 + (y - py) ** 2 + (z - pz) ** 2)) / np.sqrt(
        np.pi * np.sqrt(np.pi / (2 * exp) ** 3))


def f2(x, y, z):
    gaus1 = 0
    for i in range(len(basis['atoms'][0]['shells'][0]['con_coefficients'])):
        expi = basis['atoms'][0]['shells'][0]['p_exponents'][i]
        coefi = basis['atoms'][0]['shells'][0]['con_coefficients'][i]
        gaus1 += coefi * g1(x, y, z, -1, 0, 0, expi)

    return gaus1/np.sqrt(13.161381648762054)


def f3(x, y, z):
    gaus2 = 0
    for i in range(len(basis['atoms'][0]['shells'][0]['con_coefficients'])):
        coefi = basis['atoms'][1]['shells'][0]['con_coefficients'][i]
        expi = basis['atoms'][1]['shells'][0]['p_exponents'][i]
        gaus2 += coefi * g1(x, y, z, 1, 0, 0, expi)
    return gaus2/np.sqrt(0.015284621636345874)


#f = lambda x, y, z:   f3(x, y, z)**2

f = lambda x, y, z: (f2(x, y, z) + f3(x, y, z)) ** 2

# f = lambda x, y, z:  g1(x, y, z, 0.5, 0.5, 0.5) * g1(x, y, z, 0.5, 0.5, 0.5)
x1 = lambda y, z: -10  # Lower boundary for x
x2 = lambda y, z: 10  # Upper boundary for x
y1 = lambda z: -10  # Lower boundary for y
y2 = lambda z: 10  # Upper boundary for y
z1 = -10
z2 = 10
print(basis)
integral = integrate.tplquad(f, z1, z2, y1, y2, x1, x2)
print('integral', integral[0])

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
