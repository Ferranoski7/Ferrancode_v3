from electronic_density import ElectronicDensity
import numpy as np
from rotation.rotate import rotation

class TempFunct(object):

    def __init__(self, molecule):
        self.molecule=molecule
        self.elden=ElectronicDensity(self.molecule)
        self.atomic_value=self.elden.atomic_values
    def density_valt(self,atm_l,pos,cm=None,w=None): #gives value of the density in that point
        atomic_values=self.atomic_value

        sym_l = self.molecule.elements
        coord_a = np.array(self.molecule.coordinates)

        if cm is None:
            cm=np.array(self.elden.mass_center)

        if w is None:
            w=[1,1,1,1,1]
        # --------------------  Center of mass system ----------------------

        coord_cma=coord_a-cm

        a2au=1.889725
        coord_cm=coord_cma*a2au
        rho=0
        for i in range (len(sym_l)):
            r=(pos[0]-coord_cm[i,0])**2+(pos[1]-coord_cm[i,1])**2+(pos[2]-coord_cm[i,2])**2
            rhoi=0

            for j in range (5):
                if atomic_values[sym_l[i]][2*(j+1)+1]>0:
                    norm=(2*atomic_values[sym_l[i]][2*(j+1)]/np.pi)**(3/2)
                    rhoi=rhoi+w[j]*norm*np.exp(-2*r*atomic_values[sym_l[i]][2*(j+1)])
            rho=rho+(atomic_values[sym_l[i]][0])*rhoi

        return rho