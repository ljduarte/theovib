"""Define the "Molecule" class.
"""
import os

class Molecule:
    def __init__(self, atoms, positions, gaussian_path, aimall_path, optimised=False):
        self.atoms = atoms
        self.positions = positions
        self.gaussian_path = gaussian_path
        self.aimall_path = aimall_path

 #   def __setattr__(self, __name: str, __value: Any) -> None:
 #       pass


    def opt(self, nproc, mem,  method, base, charge, mult):
        f = open('molecule.gjf', 'w')
        print('%%nproc = %i ' %(nproc), file=f)
        print('%%mem = %i GB ' %(mem), file=f)
        print('', file=f)
        print('#opt %s %s density=current nosymm'%(method, base), file=f)
        print('', file=f)
        print('optimization', file=f)
        print('', file=f)
        print('%i %i'%(charge, mult), file=f)
        for i in range(len(self.atoms)):
            print('%s %.5f %.5f %.5f'%(self.atoms[i], self.positions[i][0], self.positions[i][1], self.positions[i][2]), file=f)
        print('', file=f)
        f.close()

        os.system(self.gaussian_path + '  molecule.gjf')




        


    
 #  def freq(nproc, mem, base, method, charge, mult):
 #   def sp(nproc, mem, base, method, charge, mult):

#    def first_derivative():

#    def second_derivative():
    




