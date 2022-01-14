from theovib.ptable import atomic_mass, constants
from theovib.aimall_tools import get_electronic
from theovib.matrices import *
from numpy import zeros, matmul
from numpy import linalg, array
"""
STEP II: Normal mode analysis 
convention: all matrix are represented by a capital letter ans small letters ate subscripts
H -> Hessian in cartesian coordinates [Hartree/Bohr]
H_iqa -> 3D matrix where each layes correspond to an IQA energy value 
M -> contain the reciprocal of the atomic mass square root   
M_t -> the transpose of M (since M is a diagonal matrix, M_t = M)
H_mw -> mass weighted Hessian [Hartree/(bohr.amu)]
L -> matrix that diagonalize H_mw. It contains theeigenvectors of H_mw
L_i -> the invese of L
L_mw -> normal coordinates = M.L
Lambda_iqa -> 3D matricex where each layer correspond to L_i multiplied by the same layer of H_iqa
Lambda -> diagonal matrix of the eigenvalues of H_mw
"""


def normal_modes(atoms, H_iqa):
    number_of_terms = len(atoms) + (len(atoms) * (len(atoms)-1))
    H = H_iqa.sum(axis=2)
    # Setting up the M matrix (construct the 1/sqrt(M) diagonal matrix) unit: [amu^-(1/2)]
    M = zeros((3*len(atoms), 3*len(atoms)))
    for i in range(3*len(atoms)):
        M[i][i] = atomic_mass[atoms[int(i/3)]]**(-0.5)
    
    #Calculate mass weighted Hessian [H_mw]: H_mw = M_t H M = M H M
    H_mw = matmul(matmul(M, H), M)
    
    #Calculate mass weighted Hessian IQA = [H_iqa_mw]
    H_iqa_mw = zeros([3*len(atoms), 3*len(atoms), number_of_terms]) #initialize the 3D matrix
    for i in range(number_of_terms):
        H_iqa_mw[:,:,i]= matmul(matmul(M, H_iqa[:,:,i]), M)
    
    #Calculate L (the eigenvectors of H_mw) and its inverse matrix
    L = linalg.eig(H_mw)[1]
    L_i = linalg.inv(L)
    
    #Use L and L_i to calculate Lambda_iqa
    Lambda_iqa = zeros([3*len(atoms), 3*len(atoms), number_of_terms]) #initialize the 3D matrix
    for i in range(number_of_terms):
        Lambda_iqa[:,:,i] =  matmul(matmul(L_i, H_iqa_mw[:,:,i]), L)
    
    L_mw = matmul(M, L) #the eigenvalues need to be ajusted by the atomic mass
    #Calculate Lambda
    Lambda = Lambda_iqa.sum(axis=2)

    # generate list of interactions keeping the reading order
    atom_label = [atoms[i] + str(i+1) for i in range(len(atoms))]
    iqa_list = ["E_intra" + "(" + i + ")" for i in atom_label]
    atom_int_label = []
    for i in range(len(atoms)):
        for j in range(i+1, len(atoms)):
            atom_int_label.append(atoms[i] + str(i+1)  + ',' + atoms[j] + str(j+1))
    iqa_list = iqa_list + ["V_cl" + '(' + i +')' for i in atom_int_label]
    iqa_list = iqa_list + ["V_xc" + '(' + i +')' for i in atom_int_label]

    #Separate the vibrational eigenvalues and calculate frequencies.
    frequencies = [] #initialize frequencies list
    frequencies_iqa = [] #initialize frequencies contributions list
    normal_coordinates =[] #initialize normal vectors list

    for i in range(len(3*atoms)):
        if Lambda[i][i] >= 0.0005: #treshold for rotacional and translational normal modes
            frequencies.append((constants['k_cm']*Lambda[i][i]/constants['k_pi'])**0.5)
            normal_coordinates.append(L_mw[:,i])
            frequencies_iqa.append(Lambda_iqa[i,i,:])
    
    return normal_coordinates, frequencies, frequencies_iqa, iqa_list

def intensities(atoms, coords, normal_coordinates, folder, delta):

    """
    STEP IIII: Eletronic properties and infrared intensities. 
    C = charge tensor
    CT = charge transfer tensor
    DP = dipolar polarization tensor
    C_prime = temporary charge derivative tensor
    X = vector of Cartesian coordinates
    X_block =square matrix of the Cartesian coordenates instances of X side by side
    A = Tranformation matrix
    AX = Hadamard product (element-wise) of the transpose of X_block and A
    """
    #initialize matrices
    C = zeros([3*len(atoms),3*len(atoms)])  #charge
    C_prime = zeros([3*len(atoms),3*len(atoms)]) #charge_derivative (temporary)
    CT = zeros([3*len(atoms),3*len(atoms)]) #charge transfer
    DP = zeros([3*len(atoms),3*len(atoms)]) #dipolar polarization

    int_files = [atoms[i].lower() + str(i+1) + ".int" for i in range(len(atoms))] #list of int files
    # generate the matrix of charge tensors
    for i in range(len(int_files)):
        charge = get_electronic(folder + "/EQ_atomicfiles/"+ int_files[i])[0]
        C[3*i,3*i] = charge
        C[3*i+1,3*i+1] = charge
        C[3*i+2,3*i+2] = charge

    # generate the CT and DP tensors matrix
    for i in range(len(3*atoms)):
        for j in range(len(atoms)):
            Electronic_A = array(get_electronic(folder + "/" + str(i) + str(i) + '_A_atomicfiles/' + int_files[j])) #Point A  
            Electronic_B = array(get_electronic(folder + "/" + str(i) + str(i) + '_B_atomicfiles/' + int_files[j])) #Ponit B
            derivative = ( Electronic_A - Electronic_B )/(2*delta)
            C_prime[3*j:3*(j+1) , i]  =  derivative[0]      
            DP[3*int(i/3), i]  =  DP[3*int(i/3), i] + derivative[1]*0.529177
            DP[3*int(i/3)+1, i] = DP[3*int(i/3)+1, i] + derivative[2]*0.529177
            DP[3*int(i/3)+2, i] = DP[3*int(i/3)+2, i] + derivative[3]*0.529177

    #write Cartesians Coordinates in the vector form
    X = zeros([3*len(atoms),1])
    n=0
    for i in range(len(atoms)):
        for j in range(3):
            X[n] = coords[i][j]
            n=n+1

    X_block = vector_to_square(X) #concatenate vector X into a square matrix
    A = gen_A_matrix(len(atoms)) #generate transformation matrx
    XA = A*X_block.transpose() # Hadamard product (element-wise)
    CT = matmul(XA, C_prime)*gen_block_identity(len(atoms)) # obtain the charge transfer tensor matrix

    Atom_tensor = C + CT+ DP # calculate the atomic tensors matrix: C + CT + DP
    intensities = [] # initialize vector to store the intensities

    #obtain intensities by multiplying the atomic tensors matrix by the normal coordinates and summing all elements
    for coordinate in normal_coordinates:
        intensities.append(constants['k_int']*(matmul(Atom_tensor, coordinate).sum())**2) #calculate infrared intensities

    return intensities, C, CT, DP



