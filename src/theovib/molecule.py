"""Define the "Molecule" class.
"""
#import os
import numpy as np


class Molecule:
    """Defines the molecule class. It stores all the relevant data for the IQA/CCTDP analysis
    """

    def __init__(self, atoms, positions):
        self.atoms = atoms
        self.positions = positions
        self.b_matrix = None
        self.hessian = None
        self.iqa_hessian = None
        self.int = None
        self.freq = None
        self.iqa_freq = None
        self.energy = None
        self.iqa_energy = None
        self.normal_coordinates = None
        self.iqa_terms = None
        self.c_tensors = None
        self.ct_tensors = None
        self.dp_tensors = None
        self.internal_hessian = None
        self.iqa_forces = None
        #self.gaussian_path = gaussian_path
        #self.aimall_path = aimall_path

#    def __setattr__(self, __name: str, __value: Any) -> None:
#       pass
    @classmethod
    def read_gaussian(cls, file):
        atoms = []
        positions = []
        f = open(file).readlines()
        for line in f:
            if len(line.split()) == 4 and all(i not in line for i in ['#', '%']):
                atoms.append(line.split()[0])
                positions.append(np.array([float(line.split()[1]), float(
                    line.split()[2]), float(line.split()[3])]))
        return(cls(atoms, positions))

#   def opt(self, nproc, mem,  method, base, charge, mult):
#       f = open('molecule.gjf', 'w')
#       print('%%nproc = %i ' %(nproc), file=f)
#       print('%%mem = %i GB ' %(mem), file=f)
#       print('', file=f)
#       print('#opt %s %s density=current nosymm'%(method, base), file=f)
#       print('', file=f)
#       print('optimization', file=f)
#       print('', file=f)
#       print('%i %i'%(charge, mult), file=f)
#       for i in range(len(self.atoms)):
#          print('%s %.5f %.5f %.5f'%(self.atoms[i], self.positions[i][0], self.positions[i][1], self.positions[i][2]), file=f)
#      print('', file=f)
#      f.close()

#      os.system(self.gaussian_path + '  molecule.gjf')

#  def freq(nproc, mem, base, method, charge, mult):
#   def sp(nproc, mem, base, method, charge, mult):

#    def first_derivative():

#    def second_derivative():
