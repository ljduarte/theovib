import numpy as np
from numpy import transpose
from theovib.ptable import *
from theovib.aimall_tools import *
import scipy.linalg

"""
Some functions to deal with matrices.
"""


def vector_to_square(vector):
    """Takes a vector and returns a diagonal matrix.

    Args:
        vector (list or array): column matrix 

    Returns:
        array: diagonal matrix
    """
    a = np.concatenate([vector, vector], axis=1)
    while a.shape[1] < vector.shape[0]:
        a = np.concatenate([a, vector], axis=1)
    return a


def gen_A_matrix(n):
    """Generates the transformation matrix A

    Args:
        n (int): number of atoms

    Returns:
        array: transformation matrix
    """
    a = np.identity(3)
    b = np.concatenate([a, a], axis=0)
    for i in range(1, n-1):
        b = np.concatenate([b, a], axis=0)
    c = np.concatenate([b, b], axis=1)
    for i in range(1, n-1):
        c = np.concatenate([c, b], axis=1)
    return c


def gen_block_identity(n):
    """Generates a block indentity matrix

    Args:
        n (int): number of atoms

    Returns:
        array: block diagonal matrix
    """
    a = np.kron(np.identity(n), np.ones([3, 3]))  # kronecker product
    return a


def invert_B(B, M):
    """Inverts the B matrix.

    Args:
        B (matrix): matrix that converts Cartesian to internal coordinates
        M (_type_): diagonal matrix of invese of atomic mass square root

    Returns:
        matrix: Inverse of B
    """
    G = np.matmul(np.matmul(B, M**2), np.transpose(B))  # Wilson's G matrix
    RHO, D = scipy.linalg.eigh(G)  # Eigenvalues and eigenvector of G

    while any(rho <= 1e-5 for rho in RHO):
        for i in range(len(RHO)):
            if RHO[i] <= 1e-5:
                print(RHO[i])
                RHO = np.delete(RHO, i)
                D = np.delete(D, i, 1)
                break
    print(len(RHO))

    RHO = np.diag(RHO)
    #G_inv =  D*RHO*D_t
    G_inv = np.matmul(np.matmul(D, np.linalg.inv(RHO)), np.transpose(D))
    B_inv = np.matmul(np.matmul(M**2, np.transpose(B)), G_inv)
    return B_inv


def hessian_from_iqa(atoms, folder, delta=0.05):
    """Generates the 3-dimensional Hessian matrix from the IQA terms.

    Args:
        atoms (list): list of atoms in the system
        delta (float): displacement value that generated the non-equilibrium 
                       geometries
        folder (list): list of _atomicfiles (aimall outputs)

    Returns:
        3D array: 3D-Hessian matrix
    """
    errors = []
    H_iqa = np.zeros([3*len(atoms), 3*len(atoms), len(atoms) +
                      (len(atoms) * (len(atoms)-1))])  # initiate the cubic matrix
    # Calculate the diagonal elements
    # Read data from the equilibrium position
    IQA_EQ = np.concatenate(get_IQA(folder + '/EQ_atomicfiles', atoms), axis=0)
    for i in range(0, 3*len(atoms)):  # loop over all XY files where X=Y
        IQA_A = np.concatenate(get_IQA(
            folder + "/" + str(i) + str(i) + '_A_atomicfiles', atoms), axis=0)  # Point A
        errors.append(get_energy_from_wfn(
            folder + "/" + str(i) + str(i) + '_A.wfn') - IQA_A.sum())
        IQA_B = np.concatenate(get_IQA(
            folder + "/" + str(i) + str(i) + '_B_atomicfiles', atoms), axis=0)  # Ponit B
        errors.append(get_energy_from_wfn(
            folder + "/" + str(i) + str(i) + '_B.wfn') - IQA_B.sum())
        derivative = (IQA_A - 2*IQA_EQ + IQA_B) / \
            (delta**2)  # Calculate the derivative
        for j in range(len(derivative)):  # loop over the IQA contributions derivatives
            # Put elements at correct position and layer, converting to the Hartree/Bohr unit.
            H_iqa[i][i][j] = derivative[j]/constants['k_bohr']

    # Calcule off-diagonal elements:
    for i in range(0, 3*len(atoms)):  # loop over all XY files where X!=Y
        for j in range(0, i):
            IQA_A = np.concatenate(get_IQA(
                folder + "/" + str(i) + str(j) + '_A_atomicfiles', atoms), axis=0)  # Point A
            errors.append(get_energy_from_wfn(
                folder + "/" + str(i) + str(j) + '_A.wfn') - IQA_A.sum())
            IQA_B = np.concatenate(get_IQA(
                folder + "/" + str(i) + str(j) + '_B_atomicfiles', atoms), axis=0)  # Point B
            errors.append(get_energy_from_wfn(
                folder + "/" + str(i) + str(j) + '_B.wfn') - IQA_B.sum())
            IQA_C = np.concatenate(get_IQA(
                folder + "/" + str(i) + str(j) + '_C_atomicfiles', atoms), axis=0)  # Point C
            errors.append(get_energy_from_wfn(
                folder + "/" + str(i) + str(j) + '_C.wfn') - IQA_C.sum())
            IQA_D = np.concatenate(get_IQA(
                folder + "/" + str(i) + str(j) + '_D_atomicfiles', atoms), axis=0)  # Point D
            errors.append(get_energy_from_wfn(
                folder + "/" + str(i) + str(j) + '_D.wfn') - IQA_D.sum())
            derivative = (IQA_A - IQA_B - IQA_C + IQA_D)/(4*delta*delta)
            # Put elements at correct position and layer, converting to the Hartree/Bohr unit.
            for l in range(len(derivative)):
                H_iqa[i][j][l] = derivative[l]/constants['k_bohr']
                H_iqa[j][i][l] = derivative[l]/constants['k_bohr']

    H = H_iqa.sum(axis=2)  # Calculate HESSIAN MATRIX
    return H, H_iqa, errors


def convert_to_internal(atoms, B, H_iqa):
    """Converts the 3-D Hessian from Cartesian to Internal coordinates

    Args:
        atoms (list): list of atoms
        B (2D array): B matrix -> converts from Cartesian to Internal coordinates
        H_iqa (3D array): 3D Hessian matrix in Cartesian coordinates 

    Returns:
        3D array: 3D Hessian in internal coordinates
    """
    M = np.zeros((3*len(atoms), 3*len(atoms)))
    for i in range(3*len(atoms)):
        M[i][i] = atomic_mass[atoms[int(i/3)]]**(-0.5)

    number_of_terms = len(atoms) + (len(atoms) * (len(atoms)-1))
    H = H_iqa.sum(axis=2)
    # Determining the G matrix
    B_inv = invert_B(B, M)
    B_inv_t = transpose(B_inv)
    H_internal = 15.5689412 * \
        np.matmul(np.matmul(B_inv_t, H), B_inv)  # B_inv_t * H * B_inv

    # Conversion factor =
    # (4.3597482 x 10^{-11} dyne-cm/Hartree) x
    # (1 Bohr/0.529177249 Ǻngstrom) x
    # (1 Bohr/0.529177249 x 10^{-8} cm) (10^3 mdyne/dyne)
    # = 15.5689412 (mdyne/Ǻ) (Hartree/Bohr2)^{-1}

    # initialize the 3D matrix
    iqa_forces = np.zeros(
        [H_internal.shape[0], H_internal.shape[1], number_of_terms])
    for i in range(number_of_terms):
        iqa_forces[:, :, i] = 15.5689412 * \
            np.matmul(np.matmul(B_inv_t, H_iqa[:, :, i]), B_inv)

    return H_internal, iqa_forces
