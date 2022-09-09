"""Define the "Molecule" class.
"""
#import os
import numpy as np

class Molecule:
    """Defines the molecule class. It stores all the relevant data for the IQA/CCTDP analysis
    """
    def __init__(self, atoms, positions, charge_mult, **kwargs):
        self.atoms = atoms
        self.positions = positions
        self.charge_mult = charge_mult 
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
#    def __setattr__(self, __name: str, __value: Any) -> None:
#       pass

    @classmethod
    def read_gaussian(cls, file):
        atoms = []
        positions = []
        f = open(file).readlines()
        for line in f:
            if len(line.split()) == 2 and all(i not in line for i in ['#', '%']):
                charge_mult = line.strip()
        for line in f:
            if len(line.split()) == 4 and all(i not in line for i in ['#', '%']):
                atoms.append(line.split()[0])
                positions.append(np.array([float(line.split()[1]), float(
                    line.split()[2]), float(line.split()[3])]))
                
        return(cls(atoms, positions, charge_mult))
    
    def gen_geometries(self, path, base_method, sigma=0.05):
        def write_geo_aa(a, atoms, coords, charge_mult, delta, base_method, output):
            R = np.zeros(3*len(atoms))
            n=0
            for i in range(len(atoms)):
                for j in range(3):
                    R[n] = coords[i][j]
                    n=n+1
            RA = np.array(R)
            RB = np.array(R)
            
            RA[a] = RA[a] + delta
            RB[a] = RB[a] - delta
            
            file = open(output + str(a) + "_" + str(a) +"_A.com", 'w')
            print('%mem = 8GB', file = file)
            print('%nproc = 8', file = file)
            print('#' + base_method + ' density=current nosym output=wfn ', file = file)
            print('', file = file)
            print(str(a) + "_" + str(a) , file = file)
            print('', file = file)
            print(charge_mult, file = file)
            for i in range(len(atoms)):
                print(atoms[i] +'    ' +str(RA[3*i]) +'\t' + str(RA[3*i+1]) +'\t' +str(RA[3*i+2]) , file = file)  
            print('', file = file)
            print(str(a) + "_" + str(a) + '_A.wfn', file = file)
            print('', file = file)
            file.close()
            
            file = open(output + str(a) + "_" + str(a) +"_B.com", 'w')
            print('%mem = 8GB', file = file)
            print('%nproc = 8', file = file)
            print('#' + base_method + ' density=current nosym output=wfn ', file = file)
            print('', file = file)
            print(str(a) + "_" + str(a) , file = file)
            print('', file = file)
            print(charge_mult, file = file)
            for i in range(len(atoms)):
                print(atoms[i] +'    ' +str(RB[3*i]) +'\t' + str(RB[3*i+1]) +'\t' +str(RB[3*i+2]) , file = file)  
            print('', file = file)
            print(str(a) + "_" + str(a) + '_B.wfn', file = file)
            print('', file = file)
            file.close()
            
            return None

        def write_geo_ab(a, b, atoms, coords, charge_mult, delta, base_method, output):
            ## Converter coordenadas num vetor
            R = np.zeros(3*len(atoms))
            n=0
            for i in range(len(atoms)):
                for j in range(3):
                    R[n] = coords[i][j]
                    n=n+1
            RA = np.array(R)
            RB = np.array(R)
            RC = np.array(R)
            RD = np.array(R)

            RA[a] = RA[a] + delta
            RA[b] = RA[b] + delta

            RB[a] = RB[a] + delta
            RB[b] = RB[b] - delta
            
            RC[a] = RC[a] - delta
            RC[b] = RC[b] + delta
            
            RD[a] = RD[a] - delta
            RD[b] = RD[b] - delta
        
            # a e b sÃ£o indices das coordenadas
            file = open(output + str(a) + "_" + str(b) +"_A.com", 'w')
            print('%mem = 8GB', file = file)
            print('%nproc = 8', file = file)
            print('#' + base_method + ' density=current nosym output=wfn ', file = file)
            print('', file = file)
            print(str(a) + "_" + str(b) , file = file)
            print('', file = file)
            print(charge_mult, file = file)
            for i in range(len(atoms)):
                print(atoms[i] +'    ' +str(RA[3*i]) +'\t' + str(RA[3*i+1]) +'\t' +str(RA[3*i+2]) , file = file)  
            print('', file = file)
            print(str(a) + "_" + str(b) + '_A.wfn', file = file)
            print('', file = file)
            file.close()
            
            file = open(output + str(a) + "_" + str(b) +"_B.com", 'w')
            print('%mem = 8GB', file = file)
            print('%nproc = 8', file = file)
            print('#' + base_method + ' density=current nosym output=wfn ', file = file)
            print('', file = file)
            print(str(a) + "_" + str(b) , file = file)
            print('', file = file)
            print(charge_mult, file = file)
            for i in range(len(atoms)):
                print(atoms[i] +'    ' +str(RB[3*i]) +'\t' + str(RB[3*i+1]) +'\t' +str(RB[3*i+2]) , file = file)  
            print('', file = file)
            print(str(a) + "_" + str(b) + '_B.wfn', file = file)
            print('', file = file)
            file.close()
            
            file = open(output + str(a) + "_" + str(b) +"_C.com", 'w')
            print('%mem = 8GB', file = file)
            print('%nproc = 8', file = file)
            print('#' + base_method + ' density=current nosym output=wfn  ', file = file)
            print('', file = file)
            print(str(a) + "_" + str(b) , file = file)
            print('', file = file)
            print(charge_mult, file = file)
            for i in range(len(atoms)):
                print(atoms[i] +'    ' +str(RC[3*i]) +'\t' + str(RC[3*i+1]) +'\t' +str(RC[3*i+2]) , file = file)  
            print('', file = file)
            print( str(a) + "_" + str(b) + '_C.wfn', file = file)
            print('', file = file)
            file.close()
            
            file = open(output + str(a) + "_" + str(b) +"_D.com", 'w')
            print('%mem = 8GB', file = file)
            print('%nproc = 8', file = file)
            print('#' + base_method + ' density=current nosym output=wfn ', file = file)
            print('', file = file)
            print(str(a) + "_" + str(b) , file = file)
            print('', file = file)
            print(charge_mult, file = file)
            for i in range(len(atoms)):
                print(atoms[i] +'    ' +str(RD[3*i]) +'\t' + str(RD[3*i+1]) +'\t' +str(RD[3*i+2]) , file = file)  
            print('', file = file)
            print(str(a) + "_" + str(b) + '_D.wfn', file = file)
            print('', file = file)
            file.close()
            
            return None
        file = open(path + "EQ.com", 'w')
        print('%mem = 8GB', file = file)
        print('%nproc = 8', file = file)
        print('#' + base_method + ' density=current nosym output=wfn', file = file)
        print('', file = file)
        print('equilibirum_geometry', file = file)
        print('', file = file)
        print(self.charge_mult, file = file)
        for i in range(len(self.atoms)):
            print(self.atoms[i] +'    ' +str(self.positions[i][0]) +'\t' +str(self.positions[i][1]) +'\t ' +str(self.positions[i][2]) , file = file)
        print('', file = file)
        print('EQ.wfn', file = file)
        print('', file = file)
        file.close()

        ### Gerar input fora da diagonal
        for a in range(1, 3*len(self.atoms), 1):
            for b in range(0,  a, 1):
                write_geo_ab(a, b, self.atoms, self.positions, self.charge_mult, sigma, base_method, path)
                
        for a in range(0, 3*len(self.atoms), 1):
            write_geo_aa(a, self.atoms, self.positions, self.charge_mult, sigma, base_method, path)
               
        return None
