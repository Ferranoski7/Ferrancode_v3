from electronic_density import ElectronicDensity
import numpy as np


class TemporalFunction(object):

    def __init__(self, molecule):
        self.molecule = molecule
        self.elden = ElectronicDensity(self.molecule)
        self.atomic_value = self.elden.atomic_values

    def density_valt(self, pos, cm=None, w=None):  # gives value of the density in that point
        atomic_values = self.atomic_value

        sym_l = self.molecule.elements
        coord_a = np.array(self.molecule.coordinates)

        if cm is None:
            cm = np.array(self.elden.mass_center())

        if w is None:
            w = [1, 1, 1, 1, 1,1,1]
        # --------------------  Center of mass system ----------------------

        coord_cma = coord_a - cm

        a2au = 1.889725
        coord_cm = coord_cma * a2au
        rho = 0
        for i in range(len(sym_l)):
            r = (pos[0] - coord_cm[i, 0]) ** 2 + (pos[1] - coord_cm[i, 1]) ** 2 + (pos[2] - coord_cm[i, 2]) ** 2
            rhoi = 0

            num_weights = self.elden.num_weights(sym_l[i])
            for j in range(num_weights):
                if atomic_values[sym_l[i]][2 * (j + 1) + 1] > 0:
                    norm = (2 * atomic_values[sym_l[i]][2 * (j + 1)] / np.pi) ** (3 / 2)
                    rhoi = rhoi + w[j] * norm * np.exp(-2 * r * atomic_values[sym_l[i]][2 * (j + 1)])
            rho = rho + (atomic_values[sym_l[i]][0]) * rhoi

        return rho

    def linear_densityt(self, init=[5, 0, 0], end=[-5, 0, 0], steps=200,
                        w=None):  # gives two outputs linear density and linear trajectory
        linear_rho = []
        length = np.sqrt((end[0] - init[0]) ** 2 + (end[1] - init[1]) ** 2 + (end[2] - init[2]) ** 2)
        trajectory = np.linspace(-length / 2, length / 2, steps)

        for i in range(steps):
            pos = [init[0] + (end[0] - init[0]) * i / (steps - 1),
                   init[1] + (end[1] - init[1]) * i / (steps - 1),
                   init[2] + (end[2] - init[2]) * i / (steps - 1)]

            linear_rho.append(self.density_valt(pos, w=w))

        linear_rho = np.array(linear_rho) * 2 * np.pi * trajectory ** 2
        return linear_rho, trajectory

    def weights(self):
        element = self.molecule.name
        atomic_values = self.elden.atomic_values
        if atomic_values[element][5] == 0:
            total = 1
        elif atomic_values[element][7] == 0:
            total = 2
        elif atomic_values[element][9] == 0:
            total = 3
        elif atomic_values[element][11] == 0:
            total = 4
        elif atomic_values[element][13] == 0:
            total = 5
        elif atomic_values[element][15] == 0:
            total = 6
        else:
            total = 7
        return total

    def gain_weights(self):
        atomic_values = self.elden.atomic_values
        element = self.molecule.name
        totalw = self.weights()
        init = [-10, 0, 0]
        end = [10, 0, 0]

        steps = 401

        spline, x = self.elden._atom_valence_spline(init, end, steps)

        valence = (atomic_values[element][0] - atomic_values[element][12]) / atomic_values[element][0]
        resta = 10000000
        itera = 10

        print(totalw)
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
                            w6 = 0
                            w7 = 0
                            w = [w1, w2, w3, w4, w5, w6, w7]
                            rho, _ = self.linear_densityt(init, end, steps, w)
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
                        w6=0
                        w7=0
                        w = [w1, w2, w3, w4, w5,w6,w7]
                        rho, _ = self.linear_densityt(init, end, steps, w)
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
                    w6 = 0
                    w7 = 0
                    w = [w1, w2, w3, w4, w5, w6, w7]
                    rho, _ = self.linear_densityt(init, end, steps, w)
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
                w6 = 0
                w7 = 0
                w = [w1, w2, w3, w4, w5, w6, w7]
                rho, _ = self.linear_densityt(init, end, steps, w)
                if np.abs(rho - spline).sum() < resta:
                    resta = np.abs(rho - spline).sum()
                    wv = w
        elif totalw == 1:
            w1 = valence
            w2 = 0
            w3 = 0
            w4 = 0
            w5 = 0
            w6 = 0
            w7 = 0
            wv = [w1, w2, w3, w4, w5, w6, w7]

        print('\n\n______________________________________________________________\n\n')
        print('The best match of Valence weights for ' + element + ' would be:')

        print('w1 = ', wv[0])
        print('w2 = ', wv[1])
        print('w3 = ', wv[2])
        print('w4 = ', wv[3])
        print('w5 = ', wv[4])
        print('w6 = ', wv[5])
        print('w7 = ', wv[6])

        return np.array(wv)
