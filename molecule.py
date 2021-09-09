import numpy as np


class Molecule(object):
    """"Creates a molecule and ta"""
    def __init__(self, elements, coordinates, name='Molecule'):
        self.elements = elements
        self.coordinates = np.array(coordinates)
        self.name=name
    def __string__(self):
        return str(self)
    def __len__(self):
        return len(self.elements)

