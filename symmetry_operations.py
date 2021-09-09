import numpy as np
from electronic_density import ElectronicDensity
from rotation.rotate import rotation


class SymOp(object):
    def __init__(self, molecule):
        self.molecule = molecule
        self.elden = ElectronicDensity(self.molecule)

    def analytic_selfsym(self, coordinates=None, axis=[0, 0, 1], angle=0, cm=None):  # Add molecule and the angle you want to check its selfsymilarity with
        atomic_values = self.elden.atomic_values
        sym_l = self.molecule.elements
        nat = len(sym_l)
        coord_a = coordinates

        if coordinates is None:
            coord_a = np.array(self.molecule.coordinates)
        if cm is None:
            cm = np.array(self.elden.mass_center())
        # --------------------  Center of mass system ----------------------
        coord_cma = coord_a - cm

        a2au = 1.889725
        coord_cm = coord_cma * a2au
        coord_cmrot = rotation(coord_cm, axis, angle)
        selfsim = 0
        for i in range(nat):
            for j in range(nat):
                for k in range(7):
                    for l in range(7):
                        a = 2 * atomic_values[sym_l[i]][2 * (k + 1)]
                        b = 2 * atomic_values[sym_l[j]][2 * (l + 1)]

                        normi = atomic_values[sym_l[i]][2 * (k + 1) + 1] * atomic_values[sym_l[i]][0] * (a / np.pi) ** (
                                3 / 2)
                        normj = atomic_values[sym_l[j]][2 * (l + 1) + 1] * atomic_values[sym_l[j]][0] * (b / np.pi) ** (
                                3 / 2)

                        coef = (np.pi / (a + b)) ** (3 / 2)
                        expn = -1.0 * a * b / (a + b)
                        r2 = (coord_cm[i, 0] - coord_cmrot[j, 0]) ** 2 + (coord_cm[i, 1] - coord_cmrot[j, 1]) ** 2 + (
                                coord_cm[i, 2] - coord_cmrot[j, 2]) ** 2
                        ex = np.exp(expn * r2)

                        selfsim = selfsim + normi * normj * coef * ex
        return selfsim

    def selfsim(self, rho1, rho2, step=0.1):  # Self similarity calculation for two arrays
        sim = (rho1 * rho2 * step * step * step).sum()
        return sim

    def analytic_invers(self, cm=None):
        ssim = self.analytic_selfsym(cm=cm)
        coord_a = np.array(self.molecule.coordinates)
        if cm is None:
            cm = np.array(self.elden.mass_center())
        # --------------------  Center of mass system ----------------------
        coord_cma = coord_a - cm
        a2au = 1.889725
        coord_cm = coord_cma * a2au
        coord_cmrot = -coord_cm

        invers = self.analytic_selfsym(coordinates=coord_cmrot)

        return invers / ssim

    def xcy(self,X,Y,axis=[0,0,1],cm=None):#2C6
        angle=X*2*np.pi/Y

        ssim=self.analytic_selfsym(axis=axis, angle=0, cm=cm)
        xcy = self.analytic_selfsym(axis=axis, angle=angle, cm=cm)/ssim

        return xcy

    def soev(self,Y,axis=[0,0,1],cm=None):
        soev=0
        for i in range(1,Y+1):
            pas=self.xcy(i,Y,axis,cm)/Y
            soev=soev+pas
        soev=1-soev
        return soev