import matplotlib.pyplot as plt
from file_ed.read import read_txt
from electronic_density import ElectronicDensity
from molecule import Molecule
from symmetry_operations import SymOp


benzene = read_txt('benzene')
benzene_sym = SymOp(benzene)
print(benzene_sym.soev(6))
#benzene_dens = ElectronicDensity(benzene)
#benzene_dens.create_cube()
























exit()

llista = ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']
element = 'Na'
atm_l = [[element, 0, 0, 0]]
wv = [1 / 11, 0, 0, 0, 0]
wc = [atomic_values[element][3] - wv[0], atomic_values[element][5] - wv[1], atomic_values[element][7] - wv[2],
      atomic_values[element][9] - wv[3], atomic_values[element][11] - wv[4]]
init = [-10, 0, 0]
end = [10, 0, 0]
length = 0
for i in range(3):
    length = length + (end[i] - init[i]) ** 2
length = np.sqrt(length)

steps = 501

rho, x = linear_density(atm_l, init, end, steps)
valence, _ = linear_densityt(atm_l, init, end, steps, wv)
core, _ = linear_densityt(atm_l, init, end, steps, wc)




plt.figure('Separation_Core_and_Valence_' + element)
plt.plot(x, rho, label='linear density')
plt.plot(x, valence, label='Valence density')
plt.plot(x, core, label='Core density')
plt.legend()
plt.xlabel('r(A)')
plt.ylabel('Radial distribution')
plt.show()














exit()

for element in range(0):
    element = llista[element]
    atm_l = [[element, 0, 0, 0]]
    totalw = weights(element)
    init = [-10, 0, 0]
    end = [10, 0, 0]
    length = 0
    for i in range(3):
        length = length + (end[i] - init[i]) ** 2
    length = np.sqrt(length)

    steps = 401

    spline, x = lin_val_density(atm_l, init, end, steps)

    valence = (atomic_values[element][0] - atomic_values[element][12]) / atomic_values[element][0]
    resta = 10000000
    itera = 10
    i, j, k, l = 0, 0, 0, 0
    if totalw == 5:
        for i in range(itera):
            total = valence
            w1 = i * valence / (itera - 1)
            total1 = total - w1
            for j in range(itera):
                w2 = j * total1 / (itera - 1)
                total2 = total1 - w2
                for k in range(itera):
                    w3 = k * total2 / (itera - 1)
                    total3 = total2 - w3
                    for l in range(itera):
                        w4 = l * total3 / (itera - 1)
                        total4 = total3 - w4
                        w5 = total4

                        w = [w1, w2, w3, w4, w5]
                        rho, _ = linear_densityt(atm_l, init, end, steps, w)
                        if np.abs(rho - spline).sum() < resta:
                            resta = np.abs(rho - spline).sum()
                            wv = w
    elif totalw == 4:
        for i in range(itera):
            total = valence
            w1 = i * valence / (itera - 1)
            total1 = total - w1
            for j in range(itera):
                w2 = j * total1 / (itera - 1)
                total2 = total1 - w2
                for k in range(itera):
                    w3 = k * total2 / (itera - 1)
                    w4 = total2 - w3
                    w5 = 0
                    w = [w1, w2, w3, w4, w5]
                    rho, _ = linear_densityt(atm_l, init, end, steps, w)
                    if np.abs(rho - spline).sum() < resta:
                        resta = np.abs(rho - spline).sum()
                        wv = w
    elif totalw == 3:
        for i in range(itera):
            total = valence
            w1 = i * valence / (itera - 1)
            total1 = total - w1
            for j in range(itera):
                w2 = j * total1 / (itera - 1)
                w3 = total1 - w2
                w4 = 0
                w5 = 0
                w = [w1, w2, w3, w4, w5]
                rho, _ = linear_densityt(atm_l, init, end, steps, w)
                if np.abs(rho - spline).sum() < resta:
                    resta = np.abs(rho - spline).sum()
                    wv = w


    elif totalw == 2:
        for i in range(itera):
            total = valence
            w1 = i * valence / (itera - 1)
            w2 = total - w1
            w3 = 0
            w4 = 0
            w5 = 0
            w = [w1, w2, w3, w4, w5]
            rho, _ = linear_densityt(atm_l, init, end, steps, w)
            if np.abs(rho - spline).sum() < resta:
                resta = np.abs(rho - spline).sum()
                wv = w
    elif totalw == 1:
        w1 = valence
        w2 = 0
        w3 = 0
        w4 = 0
        w5 = 0

    #    rho,x=linear_densityt(atm_l,init,end,steps,wv)

    #    plt.figure('Valence')
    #    plt.plot(x,rho)
    #    plt.plot(x,spline)
    #    plt.show()

    print('\n\n______________________________________________________________\n\n')
    print('The best match of Valence weights for ' + element + ' would be:')

    print('w1 = ', wv[0])
    print('w2 = ', wv[1])
    print('w3 = ', wv[2])
    print('w4 = ', wv[3])
    print('w5 = ', wv[4])


