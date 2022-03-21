import numpy as np
from numpy import transpose
from theovib.ptable import *
from theovib.aimall_tools import *
import scipy.linalg

"""
Some functions to deal with matrices.
"""
def vector_to_square(vector):
    """
    INPUT: Column matrix 
    OUTPUT: Square matrix
    """
    a = np.concatenate([vector, vector], axis=1)
    while a.shape[1] < vector.shape[0]:
        a = np.concatenate([a, vector], axis=1)
    return a


def gen_A_matrix(n):
    """
    INPUT: number of atoms 
    OUTPUT: transformation matrix
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
    """
    INPUT: number of atoms 
    OUTPUT:block diagonal matrix 
    """
    a = np.kron(np.identity(n), np.ones([3, 3]))  # kronecker product
    return a


# def invert_B(B, M):
#     """
#     G: Wilson's G matrix
#     RHO: Diagonal matrix of eigenvalues of G
#     D = eingenvectors of G
#     G_inv = D*RHO*D_t
#     B_inv = M*B_t*G_inv
#     """
#     G = np.matmul(np.matmul(B, M**2), np.transpose(B))
#     RHO, D = np.linalg.eig(G)  # eigenvalues and eigenvector of G
#     repeat = False
#     redundancies = []
#     for i in range(len(RHO)):
#         if RHO[i] <= 0.00000001:
#             print(RHO[i])
#             B = np.delete(B, i, 0)
#             redundancies.append(i)
#             repeat = True
#     if repeat == True:
#         G = np.matmul(np.matmul(B, M**2), np.transpose(B))
#         RHO, D = np.linalg.eig(G)  # eigenvalues and eigenvector of G
#         RHO = np.diag(RHO)
#         G_inv = np.matmul(np.matmul(D, np.linalg.inv(RHO)), np.transpose(D))
#         B_inv = np.matmul(np.matmul(M**2, np.transpose(B)), G_inv)
#     else:
#         RHO = np.diag(RHO)
#         G_inv = np.matmul(np.matmul(D, np.linalg.inv(RHO)), np.transpose(D))
#         B_inv = np.matmul(np.matmul(M**2, np.transpose(B)), G_inv)
#     return B_inv, redundancies

    
def invert_B(B, M):
    """
    G: Wilson's G matrix
    RHO: Diagonal matrix of eigenvalues of G
    D = eingenvectors of G
    G_inv = D*RHO*D_t
    B_inv = M*B_t*G_inv
    """
    G = np.matmul(np.matmul(B, M**2), np.transpose(B))
    RHO, D = scipy.linalg.eigh(G)  # eigenvalues and eigenvector of G
#    RHO2 = np.linalg.eigh(G)[0]
#    RHO3 = scipy.linalg.eig(G)[0]
#    RHO4 = np.linalg.eig(G)[0]

    # print('*********************************************************RHO1')
    print(RHO)
    # print('*********************************************************RHO2')
    # print(RHO2)
    # print('*********************************************************RHO3')
    # print(RHO3)  
    # print('*********************************************************RHO4')
    # print(RHO4)
    while any(rho <= 1e-8 for rho in RHO):
        for i in range(len(RHO)):
            if RHO[i] <= 1e-8:
                print(RHO[i])
                RHO = np.delete(RHO, i)
                D = np.delete(D, i, 1)
                break
    print(len(RHO))

    RHO = np.diag(RHO)
    G_inv = np.matmul(np.matmul(D, np.linalg.inv(RHO)), np.transpose(D))
    B_inv = np.matmul(np.matmul(M**2, np.transpose(B)), G_inv)
    return B_inv

def hessian_from_iqa(atoms, delta, folder):
    errors = []
    H_iqa = np.zeros([3*len(atoms), 3*len(atoms), len(atoms) + (len(atoms) * (len(atoms)-1))]) # initiate the cubic matrix
    # Calculate the diagonal elements
    IQA_EQ = np.concatenate(get_IQA(folder + '/EQ_atomicfiles', atoms), axis=0) # Read data from the equilibrium position
    for i in range(0,3*len(atoms)): #loop over all XY files where X=Y
        IQA_A = np.concatenate(get_IQA(folder + "/" + str(i) + str(i) + '_A_atomicfiles' , atoms), axis=0) #Point A
        errors.append(get_energy_from_wfn(folder + "/" + str(i) + str(i) + '_A.wfn') - IQA_A.sum())
        IQA_B = np.concatenate(get_IQA(folder + "/" + str(i) + str(i) + '_B_atomicfiles', atoms), axis=0) #Ponit B
        errors.append(get_energy_from_wfn(folder + "/" + str(i) + str(i) + '_B.wfn') - IQA_B.sum())
        derivative = (IQA_A - 2*IQA_EQ + IQA_B)/(delta**2) #Calculate the derivative    
        for j in range(len(derivative)): #loop over the IQA contributions derivatives
            H_iqa[i][i][j] = derivative[j]/constants['k_bohr'] #Put elements at correct position and layer, converting to the Hartree/Bohr unit. 

    # Calcule off-diagonal elements:    
    for i in range(0,3*len(atoms)): #loop over all XY files where X!=Y
        for j in range(0, i):     
            IQA_A = np.concatenate(get_IQA(folder + "/" + str(i) + str(j) + '_A_atomicfiles', atoms), axis=0) #Point A
            errors.append(get_energy_from_wfn(folder + "/" + str(i) + str(j) + '_A.wfn') - IQA_A.sum())
            IQA_B = np.concatenate(get_IQA(folder + "/" + str(i) + str(j) + '_B_atomicfiles', atoms), axis=0) #Point B
            errors.append(get_energy_from_wfn(folder + "/" + str(i) + str(j) + '_B.wfn') - IQA_B.sum())
            IQA_C = np.concatenate(get_IQA(folder + "/" + str(i) + str(j) + '_C_atomicfiles', atoms), axis=0) #Point C
            errors.append(get_energy_from_wfn(folder + "/" + str(i) + str(j) + '_C.wfn') - IQA_C.sum())
            IQA_D = np.concatenate(get_IQA(folder + "/" + str(i) + str(j) + '_D_atomicfiles', atoms), axis=0) #Point D
            errors.append(get_energy_from_wfn(folder + "/" + str(i) + str(j) + '_D.wfn') - IQA_D.sum())
            derivative = (IQA_A - IQA_B -IQA_C + IQA_D)/(4*delta*delta)    
            for l in range(len(derivative)): #Put elements at correct position and layer, converting to the Hartree/Bohr unit.
                H_iqa[i][j][l] = derivative[l]/constants['k_bohr']
                H_iqa[j][i][l] = derivative[l]/constants['k_bohr']

    H =  H_iqa.sum(axis=2) #Calculate HESSIAN MATRIX
    return H, H_iqa, errors

def convert_to_internal(atoms, B, H_iqa):

    """        
    STEP V: Internal coordinates         
    Internal Coordinates = B*X
    Force constants in internal coordinates:
    H_internal = B_inv_t * H * B_inv
    """        
    M = np.zeros((3*len(atoms), 3*len(atoms)))
    for i in range(3*len(atoms)):
        M[i][i] = atomic_mass[atoms[int(i/3)]]**(-0.5)

    
    number_of_terms = len(atoms) + (len(atoms) * (len(atoms)-1))
    H = H_iqa.sum(axis=2)
    #Determining the G matrix  
    B_inv = invert_B(B, M)
    B_inv_t = transpose(B_inv)
    H_internal = 15.5689412*np.matmul(np.matmul(B_inv_t, H), B_inv)
    

    """
    Conversion factor = 
    (4.3597482 x 10^{–11} dyne–cm/Hartree) x
     (1 Bohr/0.529177249 Ǻngstrom) x
     (1 Bohr/0.529177249 x 10^{–8} cm) (10^3 mdyne/dyne) 
     = 15.5689412 (mdyne/Ǻ) (Hartree/Bohr2)^{–1} 
    """
    iqa_forces = np.zeros([H_internal.shape[0],H_internal.shape[1] ,number_of_terms]) #initialize the 3D matrix
    for i in range(number_of_terms):
        iqa_forces[:,:,i]=15.5689412*np.matmul(np.matmul(B_inv_t, H_iqa[:,:,i]),B_inv)
    
    return H_internal, iqa_forces