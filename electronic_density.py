import numpy as np
from scipy.interpolate import CubicSpline


class ElectronicDensity(object):
    def __init__(self, molecule):
        self.molecule = molecule

        self.atomic_values = {
            "Titles": ['Atomic Number', 'Atomic weight', 'a1', 'w1', 'a2', 'w2', 'a3', 'w3', 'a4', 'w4', 'a5', 'w5',
                       'core electrons', 'wv1', 'wv2', 'wv3', 'wv4', 'wv5'],
            "H": [1, 1, 0.440920529854, 0.446848399957, 0.148523571410, 0.260661566642, 0.148522754226, 0.152052580541,
                  0.830533399065, 0.0886973386489, 1.55545886958, 0.0502242590213, 1, 0, 0, 0, 0, 0],
            "He": [2, 4, 0.705128316152, 0.828896903794, 4.60226167419, 0.171103096206, 1, 0,
                   1, 0, 1, 0, 2, 0, 0, 0, 0, 0],
            "Li": [3, 7, 2.42959379596, 0.496667117048, 0.0931481410958, 0.428950159018, 14.6896849265, 0.0743827239327,
                   1, 0, 1, 0, 2, 0, 1 / 3, 0, 0, 0],
            "Be": [4, 9, 4.70365863654, 0.375887797196, 0.133243850840, 0.569822506390, 27.9324424462, 0.0542896964145,
                   1, 0, 1, 0, 2, 0, 1 / 2, 0, 0, 0],
            "B": [5, 11, 7.69672609374, 0.300184331682, 0.209275915739, 0.657081391054, 45.2385931510, 0.0427342772638,
                  1, 0, 1, 0, 2, 0, 0.6, 0, 0, 0],
            "C": [6, 12, 0.356462048751, 0.481662583478, 6.74589838364, 0.165621922856, 20.0053616858, 0.126630091695,
                  97.0508229097, 0.0169326970446, 0.122734548099, 0.209152704926, 2, 4 / 9, 0, 0, 0, 2 / 9],
            "N": [7, 14, 0.558690868187, 0.428044393139, 26.3091078048, 0.119405132277, 132.716385327, 0.0149293216485,
                  0.213410401209, 0.306224338493, 8.67117752141, 0.131396814443, 2, 0.396825, 0, 0, 0.31746, 0],
            "O": [8, 16, 0.764949764007, 0.462949674244, 36.3980482778, 0.0972634518752, 178.279439955, 0.0124323059053,
                  0.258345655165, 0.306747136544, 12.1009178518, 0.120607431432, 2, 0.416666, 0, 0, 1 / 3, 0],
            "F": [9, 19, 0.953231108644, 0.521113130773, 0.288290993131, 0.272227674853, 14.3826569656, 0.102468771030,
                  44.6764228163, 0.0927110692516, 223.415128313, 0.0114793540601, 2, 0, 7 / 9, 0, 0, 0],
            "Ne": [10, 20, 0.954374750854E+00, 0.828487153815E+00, 0.182218394930E+03, 0.229002856433E-01,
                   0.320674536956E+02, 0.148612560541E+00,
                   1, 0, 1, 0, 2, 0, 0.8, 0, 0, 0],
            "Na": [11, 23, 0.108854352851E+01, 0.825554838994E+00, 0.242558589286E+02, 0.112378067204E+00,
                   0.933932018858E+02, 0.577146198741E-01,
                   0.529038388855E+03, 0.435247392718E-02, 1, 0, 10, 0, 0, 0, 0, 0],
            "Mg": [12, 24, 0.218569903740E+01, 0.520863441795E+00, 0.301497473367E+00, 0.342512252791E+00,
                   0.490777626016E+02, 0.117633186674E+00,
                   0.267399143111E+03, 0.189911187397E-01, 1, 0, 10, 0, 0.166666, 0, 0, 0],
            "Al": [13, 27, 0.261388211922E+01, 0.523097525672E+00, 0.214149737434E+00, 0.351056420809E+00,
                   0.588955805607E+02, 0.109087878782E+00,
                   0.323698543235E+03, 0.167581747358E-01, 1, 0, 10, 0, 0.2307, 0, 0, 0],
            "Si": [14, 28, 0.316910421369E+01, .499707795463E+00, 0.204138192913E+00, 0.384254907417E+00,
                   0.694191161972E+02, 0.100765099351E+00,
                   0.381267247680E+03, 0.152721977687E-01, 1, 0, 10, 0, 0.2857, 0, 0, 0],
            "P": [15, 31, 0.361451144211E+01, 0.478831971088E+00, 0.203232522123E+00, 0.407091911330E+00,
                  0.642302360554E+02, 0.852748018285E-01,
                  0.241672185606E+03, 0.270948640067E-01, 0.140780726438E+04, 0.170645174514E-02, 10, 0, 1 / 3, 0, 0,
                  0],
            "S": [16, 32, 0.432801783704E+01, 0.450261917524E+00, 0.234919053755E+00, 0.443771174997E+00,
                  0.748183927869E+02, 0.797906413585E-01,
                  0.281418786843E+03, 0.246258047332E-01, 0.163076922864E+04, 0.155046138826E-02, 10, 0, 0.3375, 0, 0,
                  0.0375],
            "Cl": [17, 35, 0.512735107009E+01, 0.423533567847E+00, 0.274578692637E+00, 0.477496012097E+00,
                   0.856745844595E+02, 0.745409754409E-01, 0.319382224252E+03, 0.229544216482E-01, 1, 0, 10, 0, 0.4117,
                   0, 0, 0],
            "Ar": [18, 40, 0.630553641999E+01, 0.395500421392E+00, 0.334230771034E+00, 0.516912377231E+00,
                   0.121409292169E+03, 0.764819814552E-01,
                   0.663530875228E+03, 0.111052199211E-01, 1, 0, 10, 0, 4 / 9, 0, 0, 0],

        }

    def _weights(self, element):

        atomic_values = self.atomic_values
        weights = np.array([atomic_values[element][3], atomic_values[element][5], atomic_values[element][7],
                                 atomic_values[element][9], atomic_values[element][11]])
        return weights

    def _valence_weights(self,element):
        atomic_values=self.atomic_values
        valence_weights=np.array([atomic_values[element][13],atomic_values[element][14],atomic_values[element][15],atomic_values[element][16],atomic_values[element][17]])

        return valence_weights

    def _core_weights(self, element):
        weights=self._weights(element)
        valence_weights= self._valence_weights(element)
        core_weights=weights-valence_weights
        return core_weights

    def mass_center(self):
        atomic_values=self.atomic_values
        sym_l = np.array(self.molecule.elements)
        coord_a = np.array(self.molecule.coordinates)


        sym_m = []
        for i in range(len(sym_l)):
            sym_m.append(atomic_values[sym_l[i]][1])

        sym_m=np.array(sym_m)
        cm=[np.sum(coord_a[:,0]*sym_m[:])/np.sum(sym_m), np.sum(coord_a[:,1]*sym_m[:])/np.sum(sym_m), np.sum(coord_a[:,2]*sym_m[:])/np.sum(sym_m)]

        return cm

    def density_val(self, pos, center=None):
        atomic_values = self.atomic_values
        sym_l = self.molecule.elements
        coord_a = np.array(self.molecule.coordinates)

        sym_m = []
        sym_num = []

        for i in range(len(sym_l)):
            sym_m.append(atomic_values[sym_l[i]][1])
            sym_num.append(atomic_values[sym_l[i]][0])
        if center is None:
            center = self.mass_center()

        coord_cma = coord_a - center

        a2au = 1.889725
        coord_cm = coord_cma * a2au
        rho = 0
        for i in range(len(sym_l)):
            r = (pos[0] - coord_cm[i, 0]) ** 2 + (pos[1] - coord_cm[i, 1]) ** 2 + (pos[2] - coord_cm[i, 2]) ** 2
            rhoi = 0

            for j in range(5):
                if atomic_values[sym_l[i]][2 * (j + 1) + 1] > 0:
                    norm = (2 * atomic_values[sym_l[i]][2 * (j + 1)] / np.pi) ** (3 / 2)
                    rhoi = rhoi + norm * atomic_values[sym_l[i]][2 * (j + 1) + 1] * np.exp(
                        -2 * r * atomic_values[sym_l[i]][2 * (j + 1)])
            rho = rho + (atomic_values[sym_l[i]][0]) * rhoi

        return rho

    def linear_density(self, init=[-5, 0, 0], end=[5, 0, 0], steps=200):
        atomic_values = self.atomic_values
        lrho = []
        leng = np.sqrt((end[0] - init[0]) ** 2 + (end[1] - init[1]) ** 2 + (end[2] - init[2]) ** 2)
        trajectory = np.linspace(-leng / 2, leng / 2, steps)

        for i in range(steps):
            pos = [init[0] + (end[0] - init[0]) * i / (steps - 1),
                   init[1] + (end[1] - init[1]) * i / (steps - 1),
                   init[2] + (end[2] - init[2]) * i / (steps - 1)]

            lrho.append(self.density_val(pos))

        lrho = np.array(lrho) * 2 * np.pi * trajectory ** 2
        return lrho, trajectory

    def valence_dens(self,pos,center=None):
        atomic_values = self.atomic_values
        sym_l = self.molecule.elements
        coord_a = np.array(self.molecule.coordinates)

        if center is None:
            center = self.mass_center()



        coord_cma = coord_a - center

        a2au = 1.889725
        coord_cm = coord_cma * a2au
        rho = 0
        for i in range(len(sym_l)):
            r = (pos[0] - coord_cm[i, 0]) ** 2 + (pos[1] - coord_cm[i, 1]) ** 2 + (pos[2] - coord_cm[i, 2]) ** 2
            rhoi = 0

            for j in range(5):
                if atomic_values[sym_l[i]][2 * (j + 1) + 1] > 0:
                    norm = (2 * atomic_values[sym_l[i]][2 * (j + 1)] / np.pi) ** (3 / 2)
                    rhoi = rhoi + norm * atomic_values[sym_l[i]][j+13] * np.exp(
                        -2 * r * atomic_values[sym_l[i]][2 * (j + 1)])
            rho = rho + (atomic_values[sym_l[i]][0]) * rhoi

        return rho

    def linear_valence_density(self, init=[-5, 0, 0], end=[5, 0, 0], steps=200):
        atomic_values = self.atomic_values
        lrho = []
        leng = np.sqrt((end[0] - init[0]) ** 2 + (end[1] - init[1]) ** 2 + (end[2] - init[2]) ** 2)
        trajectory = np.linspace(-leng / 2, leng / 2, steps)

        for i in range(steps):
            pos = [init[0] + (end[0] - init[0]) * i / (steps - 1),
                   init[1] + (end[1] - init[1]) * i / (steps - 1),
                   init[2] + (end[2] - init[2]) * i / (steps - 1)]

            lrho.append(self.valence_dens(pos))

        lrho = np.array(lrho) * 2 * np.pi * trajectory ** 2
        return lrho, trajectory

    def _atom_valence_spline(self, init=[-5, 0, 0], end=[5, 0, 0], steps=401):
        atomic_values = self.atomic_values
        y2, x = self.linear_density(init, end, steps)
        atm_l = self.molecule.elements
        element=str(atm_l[0])
        length = 0
        for i in range(3):
            length = length + (end[i] - init[i]) ** 2
        length = np.sqrt(length)
        dery2 = []
        dery2.append(0)
        for i in range(1, len(y2)):
            dery2.append((y2[i] - y2[i - 1]) / (x[i] - x[i - 1]))
        y3 = []
        x2 = []
        cont = 0
        counter = 0
        coef = 1
        mini = np.min(y2)

        for i in range(1, len(x) - 1):
            if cont == 0:
                if (y2[i - 1] > y2[i]) and (y2[i] < y2[i + 1]) and y2[i] > mini:
                    coef = 1 ^ coef
                    ymin = y2[i]
                    xmin = x[i]
                    cont = cont + 1
                elif (dery2[i - 1] > dery2[i]) and (dery2[i] < dery2[i + 1]) and y2[i] > mini:
                    coef = 1 ^ coef
                    ymin = y2[i]
                    xmin = x[i]
                    cont = cont + 2
            elif cont == 1:
                if (y2[i - 1] > y2[i]) and (y2[i] < y2[i + 1]) and y2[i] > mini and y2[i] / ymin < 1.2 and x[i] != xmin:
                    coef = 1 ^ coef
            elif cont == 2:
                if (dery2[i - 1] < dery2[i]) and (dery2[i] > dery2[i + 1]) and y2[i] > mini and y2[i] / ymin < 1.2 and \
                        x[i] != xmin:
                    coef = 1 ^ coef

            y3.append(y2[i] * coef)
            x2.append(x[i] * coef)
            counter = counter + (1 ^ coef)

        for i in range(counter):
            y3.remove(0)
            x2.remove(0)

        counter = 0
        for i in range(1, len(y3) - 1):
            step = x2[2] - x2[1]
            if (x2[i + 1] - x[i]) > (2 * step) and (counter != 1):
                y3.insert(i + 1, 0)
                x2.insert(i + 1, 0)
                counter = 1
        cs = CubicSpline(x2, y3)

        valence = (atomic_values[element][0] - atomic_values[element][12]) / ((cs(x)).sum() * length / steps)

        y = cs(x) * valence
        return y, x

    def grid(self,step=0.1):
        atomic_values = self.atomic_values
        sym_l = self.molecule.elements
        coord_a = np.array(self.molecule.coordinates)
        cm = self.mass_center()
        # --------------------  Center of mass system ----------------------
        coord_cma = coord_a - cm
        a2au = 1.889725
        coord_cm = coord_cma * a2au
        # ---------------------- Electronic density --------------------------

        rm = np.max(np.abs(coord_cm)) + 50 * step

        num = int(2 * rm / step) + 1

        # we create a grid that we'll use to create the electronic density centered in the center of Mass
        x_ = np.linspace(-rm, rm, num)
        y_ = np.linspace(-rm, rm, num)
        z_ = np.linspace(-rm, rm, num)

        x, y, z = np.meshgrid(x_, y_, z_, indexing='ij')

        # We'll create the electronic density for the molecule
        rho = np.zeros((x.shape[0], x.shape[1], x.shape[0]))
        for i in range(len(sym_l)):
            rmesh = (x - coord_cm[i, 0]) ** 2 + (y - coord_cm[i, 1]) ** 2 + (z - coord_cm[i, 2]) ** 2
            rhoi = 0 * rmesh

            for j in range(5):
                if atomic_values[sym_l[i]][2 * (j + 1) + 1] > 0:
                    norm = (2 * atomic_values[sym_l[i]][2 * (j + 1)] / np.pi) ** (3 / 2)
                    rhoi = rhoi + norm * atomic_values[sym_l[i]][2 * (j + 1) + 1] * np.exp(
                        -2 * rmesh * atomic_values[sym_l[i]][2 * (j + 1)])
            rho = rho + (atomic_values[sym_l[i]][0]) * rhoi
        return rho

    def create_cube(self,step=0.1):
        sysl=self.molecule.name
        atomic_values = self.atomic_values
        sym_l = np.array(self.molecule.elements)
        coord_a = np.array(self.molecule.coordinates)
        # ---------------------- Centre de masses i tal --------------------------
        cm=self.mass_center()
        # --------------------  Center of mass system ----------------------

        coord_cma = coord_a - cm

        # ---------------------- Electronic density --------------------------
        a2au = 1.889725
        coord_cm = coord_cma * a2au

        rm = np.max(np.abs(coord_cm)) + 50 * step

        num = int(2 * rm / step) + 1
        rho = self.grid(step)
        print(np.shape(rho))
        # ---------------------- Creating the .cube file -------------------------
        f = open(sysl + ".cube", "w")
        f.write("Header \n OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z \n")
        # the center of the plot will be displaced by rm in the .cube file for the volumetric data
        cen = rm

        f.write(str(len(sym_l)) + '\t' + str(0) + '\t' + str(0) + '\t' + str(0) + '\n')

        # We define the plot size and unitary vectors used (ortogonal in this case)
        f.write(str(num) + '\t' + str(step) + '\t0\t0\n')
        f.write(str(num) + '\t0\t' + str(step) + '\t0\n')
        f.write(str(num) + '\t0\t0\t' + str(step) + '\n')

        # Placing of the different atoms
        for i in range(len(sym_l)):
            f.write(str(atomic_values[sym_l[i]][0]) + '\t0\t' + str(coord_cm[i, 0] + cen) + '\t' + str(
                coord_cm[i, 1] + cen) + '\t' + str(coord_cm[i, 2] + cen) + '\n')

        # Applying the .cube the volumetric data format to the electronic density
        for i in range(num):
            for j in range(num):
                for k in range(num):
                    f.write(str(rho[i, j, k]) + '\t')
                    if ((k) % 6) == 5:
                        f.write('\n')

                f.write('\n')

        f.close()
        return
