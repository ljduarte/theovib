from theovib.ptable import atomic_mass, constants
from theovib.aimall_tools import get_electronic
from theovib.matrices import *
from numpy import zeros, matmul
from numpy import linalg, array


def normal_modes(atoms, H_iqa):
    """Calculates the normal modes of vibration from the IQA 3D-Hessian matrix

    Args:
        atoms (list): list of atoms in the system
        H_iqa (3D array): 3D Hessian calculated with IQA contributions

    Returns:
        2D array: normal modes matrix. Each column contrains the displacement of atoms for a normal coordinate
        list: frequencies of vibration
        3D array: IQA partitioning of frequencies
        list: list of IQA terms
    """
    number_of_terms = len(atoms) + (len(atoms) * (len(atoms)-1))
    H = H_iqa.sum(axis=2)
    # Setting up the M matrix (construct the 1/sqrt(M) diagonal matrix) unit: [amu^-(1/2)]
    M = zeros((3*len(atoms), 3*len(atoms)))
    for i in range(3*len(atoms)):
        M[i][i] = atomic_mass[atoms[int(i/3)]]**(-0.5)

    # Calculate mass weighted Hessian [H_mw]: H_mw = M_t H M = M H M
    H_mw = matmul(matmul(M, H), M)

    # Calculate mass weighted Hessian IQA = [H_iqa_mw]
    # initialize the 3D matrix
    H_iqa_mw = zeros([3*len(atoms), 3*len(atoms), number_of_terms])
    for i in range(number_of_terms):
        H_iqa_mw[:, :, i] = matmul(matmul(M, H_iqa[:, :, i]), M)

    # Calculate L (the eigenvectors of H_mw) and its inverse matrix
    L = linalg.eig(H_mw)[1]
    L_i = linalg.inv(L)

    # Use L and L_i to calculate Lambda_iqa
    # initialize the 3D matrix
    Lambda_iqa = zeros([3*len(atoms), 3*len(atoms), number_of_terms])
    for i in range(number_of_terms):
        Lambda_iqa[:, :, i] = matmul(matmul(L_i, H_iqa_mw[:, :, i]), L)

    # the eigenvalues need to be ajusted by the atomic mass
    L_mw = matmul(M, L)
    # Calculate Lambda
    Lambda = Lambda_iqa.sum(axis=2)

    # generate list of interactions keeping the reading order
    atom_label = [atoms[i] + str(i+1) for i in range(len(atoms))]
    iqa_list = ["E_intra" + "(" + i + ")" for i in atom_label]
    atom_int_label = []
    for i in range(len(atoms)):
        for j in range(i+1, len(atoms)):
            atom_int_label.append(
                atoms[i] + str(i+1) + ',' + atoms[j] + str(j+1))
    iqa_list = iqa_list + ["V_cl" + '(' + i + ')' for i in atom_int_label]
    iqa_list = iqa_list + ["V_xc" + '(' + i + ')' for i in atom_int_label]

    # Separates the vibrational eigenvalues and calculate frequencies.
    frequencies = []  # initializes frequencies list
    frequencies_iqa = []  # initializes frequencies contributions list
    normal_coordinates = []  # initializes normal vectors list

    for i in range(len(3*atoms)):
        if Lambda[i][i] >= 0.0005:  # treshold for rotational and translational normal modes
            frequencies.append(
                (constants['k_cm']*Lambda[i][i]/constants['k_pi'])**0.5)
            normal_coordinates.append(L_mw[:, i])
            frequencies_iqa.append(Lambda_iqa[i, i, :])

    return normal_coordinates, frequencies, frequencies_iqa, iqa_list


def intensities(atoms, coords, normal_coordinates, folder, delta):
    """Calculates the infrared intensities and CCTDP contributions from AIM charges and dipoles 

    Args:
        atoms (list): list of .int files (AIMAll outputs)
        coords (2D array): Cartesian coordinates
        normal_coordinates (2D array): normal coordinates matrix
        folder (list): list of AIMAll _atomicfiles folders
        delta (float): displacement value that generated the non-equilibrium geometries

    Returns:
        list: infrared intensities
        2D array: Charge (C) contribution
        2D array: Charge-Transfer (CT) contribution
        2D array: Dipolar Polarization (DP) contribution
    """
    # initializes matrices
    C = zeros([3*len(atoms), 3*len(atoms)])  # charge
    # charge_derivative (temporary)
    C_prime = zeros([3*len(atoms), 3*len(atoms)])
    CT = zeros([3*len(atoms), 3*len(atoms)])  # charge transfer
    DP = zeros([3*len(atoms), 3*len(atoms)])  # dipolar polarization

    # list of int files
    int_files = [atoms[i].lower() + str(i+1) +
                 ".int" for i in range(len(atoms))]
    # generate the matrix of charge tensors
    for i in range(len(int_files)):
        charge = get_electronic(folder + "/EQ_atomicfiles/" + int_files[i])[0]
        C[3*i, 3*i] = charge
        C[3*i+1, 3*i+1] = charge
        C[3*i+2, 3*i+2] = charge

    # generate the CT and DP tensors matrix
    for i in range(len(3*atoms)):
        for j in range(len(atoms)):
            Electronic_A = array(get_electronic(
                folder + "/" + str(i) + '_' + str(i) + '_A_atomicfiles/' + int_files[j]))  # Point A
            Electronic_B = array(get_electronic(
                folder + "/" + str(i) + '_' +  str(i) + '_B_atomicfiles/' + int_files[j]))  # Ponit B
            derivative = (Electronic_A - Electronic_B)/(2*delta)
            C_prime[3*j:3*(j+1), i] = derivative[0]
            DP[3*int(i/3), i] = DP[3*int(i/3), i] + derivative[1]*0.529177
            DP[3*int(i/3)+1, i] = DP[3*int(i/3)+1, i] + derivative[2]*0.529177
            DP[3*int(i/3)+2, i] = DP[3*int(i/3)+2, i] + derivative[3]*0.529177

    # write Cartesians Coordinates in the vector form
    X = zeros([3*len(atoms), 1])
    n = 0
    for i in range(len(atoms)):
        for j in range(3):
            X[n] = coords[i][j]
            n = n+1

    X_block = vector_to_square(X)  # concatenate vector X into a square matrix
    A = gen_A_matrix(len(atoms))  # generate transformation matrx
    XA = A*X_block.transpose()  # Hadamard product (element-wise)
    # obtain the charge transfer tensor matrix
    CT = matmul(XA, C_prime)*gen_block_identity(len(atoms))

    Atom_tensor = C + CT + DP  # calculate the atomic tensors matrix: C + CT + DP
    D = gen_D_matrix(len(atoms))
    intensities = []  # initialize vector to store the intensities

    # obtain intensities by multiplying the atomic tensors matrix by the normal coordinates and summing all elements
    for coordinate in normal_coordinates:
        # calculate infrared intensities
        PQ=matmul(D, matmul(Atom_tensor, coordinate))
        intensities.append(
            constants['k_int']*matmul(transpose(PQ), PQ))

    return intensities, C, CT, DP
