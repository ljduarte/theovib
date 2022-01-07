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
    return


def invert_B(B, M):
    """
    G: Wilson's G matrix
    RHO: Diagonal matrix of eigenvalues of G
    D = eingenvectors of G
    G_inv = D*RHO*D_t
    B_inv = M*B_t*G_inv
    """
    G = np.matmul(np.matmul(B, M), np.transpose(B))
    RHO, D = np.linalg.eig(G)  # eigenvalues and eigenvector of G
    repeat = False
    for i in range(len(RHO)):
        if RHO[i] < 0.1:
            del B[i]
            print('Redundant coordinate:', i)
            repeat = True
    if repeat == True:
        G = np.matmul(np.matmul(B, M), np.transpose(B))
        RHO, D = np.linalg.eig(G)  # eigenvalues and eigenvector of G
        RHO = np.diag(RHO)
        G_inv = np.matmul(np.matmul(D, np.linalg.inv(RHO)), np.transpose(D))
        B_inv = np.matmul(np.matmul(M, np.transpose(B)), G_inv)
    else:
        RHO = np.diag(RHO)
        G_inv = np.matmul(np.matmul(D, np.linalg.inv(RHO)), np.transpose(D))
        B_inv = np.matmul(np.matmul(M, np.transpose(B)), G_inv)
    return B_inv


