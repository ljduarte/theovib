if __name__ == '__main__':
    
    from theovib.molecule import *
    from theovib.internal import *
    from theovib.matrices import *
    from theovib.ir import *
    from theovib.output import *
    import sys 
    import pandas as pd
    
    input_data = Input.read_text('input.txt')
    molecule = Molecule.read_gaussian(input_data.folder + '/EQ.com')
    molecule.energy = get_energy_from_wfn(input_data.folder +'/EQ.wfn')
    molecule.iqa_energy = get_IQA(input_data.folder +'/EQ_atomicfiles', molecule.atoms)
    b_matrix = []
    for coord in input_data.bond:
        b_matrix.append(bond(molecule.positions, coord[0], coord[1]))
    for coord in input_data.angle:
        b_matrix.append(angle(molecule.positions, coord[0], coord[1], coord[2]))    
    for coord in input_data.torsion:
        b_matrix.append(torsion(molecule.positions, coord[0], coord[1], coord[2], coord[3]))
    for coord in input_data.linear_angle:
        b_matrix.append(linear(molecule.positions, coord[0], coord[1]))
    for coord in input_data.oop_wag:
        b_matrix.append(wag(molecule.positions, coord[0], coord[1], coord[2], coord[3]))
        
    molecule.b_matrix = np.array(b_matrix)
    molecule.hessian, molecule.iqa_hessian, errors = hessian_from_iqa(molecule.atoms, input_data.delta, input_data.folder)
    molecule.normal_coordinates, molecule.freq, molecule.iqa_freq, molecule.iqa_terms = normal_modes(molecule.atoms, molecule.iqa_hessian)
    molecule.int, molecule.c_tensors, molecule.ct_tensors, molecule.dp_tensors = intensities(molecule.atoms, molecule.positions, molecule.normal_coordinates, input_data.folder, input_data.delta)
    molecule.internal_hessian, molecule.iqa_forces = convert_to_internal(molecule.atoms, molecule.b_matrix, molecule.iqa_hessian)
    
    
    original_stdout = sys.stdout
    with open('output.txt', 'w') as f:
        sys.stdout = f # Change the standard output to the file we created.
        print('########################################\n             Theovib output\n########################################\n\nA Library to Compute Infrared Properties\nWritten by L. J. Duarte\n----------------------------------------\n')
        print('MOLECULE:\t' + input_data.molecule)
        print('FOLDER:\t' + input_data.folder)
        print('')
        print('Delta:\t' + str(input_data.delta) + '  (displacement for numerical derivatives)')
        print('')
        print('#################\n Input Geometry\n#################\n')
        print('ATOMS:\t' + '\t'.join(molecule.atoms))
        print('')
        print('CARTESIAN COORDINATES:')
        print('\n'.join(['\t'.join([str("{:10.6f}".format(cell)) for cell in row]) for row in molecule.positions]))
        print('')
        print('wfn ENERGY (Hartress):\t' + str(molecule.energy))
        print('IQA ENERGY (Hartress):\t' + str(sum(sum(molecule.iqa_energy))))
        print('\nHESSIAN:')
        print('\n'.join(['\t'.join([str("{:10.6f}".format(cell)) for cell in row]) for row in molecule.hessian]))
        print('')
        print('#################\n  Normal modes  \n#################\n')
        print('FREQUECIES [cm^{-1}]     :\t' + '\t'.join([str("{:10.2f}".format(freq)) for freq in molecule.freq]))
        print('INTENSITIES [km mol^{-1}]:\t' + '\t'.join([str("{:10.2f}".format(inte)) for inte in molecule.int]))
        print('L MATRICES:')
        print('\n'.join(['\t'.join([str("{:10.6f}".format(cell)) for cell in row]) for row in np.transpose(molecule.normal_coordinates)]))
        print('')
        print('#################\n      CCTDP   \n#################\n')
        print('Charge tensor:')
        print('\n'.join(['\t'.join([str("{:10.6f}".format(cell)) for cell in row]) for row in molecule.c_tensors]))
        print('')
        print('Charge Transfer tensor:')
        print('\n'.join(['\t'.join([str("{:10.6f}".format(cell)) for cell in row]) for row in molecule.ct_tensors]))
        print('')
        print('Dipolar Polarization tensor:')
        print('\n'.join(['\t'.join([str("{:10.6f}".format(cell)) for cell in row]) for row in molecule.dp_tensors]))
        print('')
        print('Atomic polar tensors:')
        print('\n'.join(['\t'.join([str("{:10.6f}".format(cell)) for cell in row]) for row in molecule.c_tensors + molecule.ct_tensors + molecule.dp_tensors]))
        print('')
        print('#################\n      IQA   \n#################\n')
        print('')
        print('FREQUENCY PARTITIONING\n')
        
        for i in range(len(molecule.freq)):
            print(str("{:10.2f}".format(molecule.freq[i])))
            print('\t'.join(molecule.iqa_terms))
            print('\t'.join([str("{:10.6f}".format(cell)) for cell in molecule.iqa_freq[i]])) 
            print('')
            
        print('INTERNAL COORDINATES\n')
        print('Wilson B Matrix:')
        for coord in input_data.bond:
            print('BOND ' + str(coord) )
        for coord in input_data.angle:
            print('ANGLE ' + str(coord))    
        for coord in input_data.torsion:
            print('TORSION ' + str(coord))   
        for coord in input_data.linear_angle:
            print('ANGLE ' + str(coord))   
        for coord in input_data.oop_wag:
            print('OOP WAG ' + str(coord))
        print('')
        print('\n'.join(['\t'.join([str("{:10.6f}".format(cell)) for cell in row]) for row in molecule.b_matrix]))
        print('')
        print('Internal Hessian:')
        print('\n'.join(['\t'.join([str("{:10.6f}".format(cell)) for cell in row]) for row in molecule.internal_hessian]))
        print('')
        print('IQA Forces:')
        print('')
        for i in range(len(molecule.internal_hessian)):
            for j in range(len(molecule.internal_hessian)):
                print(str(i)+str(j))
                print('\t'.join(molecule.iqa_terms))
                print('\t'.join([str("{:10.6f}".format(cell)) for cell in molecule.iqa_forces[i][j]])) 
                print('')
        print("RMSE is: %f kJ/mol"%(2625.5*np.sqrt(sum(np.array(errors)**2)/len(errors))))
        print("MAE is: %f kJ/mol"%(2625.5*(sum(np.abs(np.array(errors)))/len(errors))))
    sys.stdout = original_stdout # Reset the standard output to its original value