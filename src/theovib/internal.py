import numpy as np
"""Functions that calculate the internal coordinates from Cartesian 
coordinates
"""
def bond(geo, a, b):
    """Calculates the bond stretch internal coordinate between a and b 
    in the molecule

    Args:
        geo (2D array): molecule geometry in Cartesian Coordinates
        a (int): label of atom A
        b (int): label of atom B

    Returns:
        1D array: a row of the B matrix
    """    
    v = (geo[a-1]-geo[b-1])/np.linalg.norm(geo[a-1]-geo[b-1])
    b_row = np.zeros(3*len(geo))
    b_row[3*a-3] = v[0]
    b_row[3*a-2] = v[1]
    b_row[3*a-1] = v[2]
    b_row[3*b-3] = -v[0]
    b_row[3*b-2] = -v[1]
    b_row[3*b-1] = -v[2]

    return b_row

def angle(geo, a, b, c):
    """Calculates the angle bending internal coordinate between atoms 
    a-b-c

    Args:
        geo (2D array): molecular geometry in Cartesian Coordinates
        a (int): label of atom A
        b (int): label of central atom B
        c (int): label of atom C

    Returns:
        1D array: a row of the B matrix
    """
    v1 = geo[a-1]-geo[b-1]
    r1 = np.linalg.norm(v1)
    v1 = v1/r1
    v2 = geo[c-1]-geo[b-1]
    r2 = np.linalg.norm(v2)
    v2 = v2/r2
    normal = np.cross(v1, v2)/np.linalg.norm(np.cross(v1, v2))
    u1 = np.cross(v1, normal)
    u1 = u1/np.linalg.norm(u1)
    u1 = u1/r1
    u2 = np.cross(normal, v2)
    u2 = u2/np.linalg.norm(u2)
    u2 = u2/r2
    b_row = np.zeros(3*len(geo))
    b_row[3*a-3] = u1[0]
    b_row[3*a-2] = u1[1]
    b_row[3*a-1] = u1[2]
    b_row[3*c-3] = u2[0]
    b_row[3*c-2] = u2[1]
    b_row[3*c-1] = u2[2]
    b_row[3*b-3] = -(u1[0]+u2[0])
    b_row[3*b-2] = -(u1[1]+u2[1])
    b_row[3*b-1] = -(u1[2]+u2[2])
 
    return b_row

def torsion(geo, a, b, c, d):
    """Calculates the torsion internal coordinate between n atoms A, atom B, 
    atom C and m atoms D

    Args:
        geo (2D array): _description_
        a (list): list of labels for the n A-type atoms
        b (int): label of atom B
        c (int): label of atom C
        d (list): list of labels for the m D-type atoms

    Returns:
        1D array: a row of the B matrix
    """
    n = len(a)
    m = len(d)
    b_row = np.zeros(3*len(geo))
    bc = geo[c-1]-geo[b-1]
    r_bc = np.linalg.norm(bc)
    bc = bc/r_bc
    cb = -bc
    vb = [0, 0, 0]
    vc = [0, 0, 0]

    for a in a:
        v1 = geo[b-1]-geo[a-1]
        r1 = np.linalg.norm(v1)
        v1 = v1/r1
        ua = np.cross(v1, bc)/(r1*n*np.linalg.norm(np.cross(v1, bc))**2)
        b_row[3*a-3] = -ua[0]
        b_row[3*a-2] = -ua[1]
        b_row[3*a-1] = -ua[2]
        vb = vb + (r_bc-r1*(np.dot(v1, bc)))*ua/r_bc
        vc = vc + np.dot(bc, v1)*r1*ua/(r_bc)

    for d in d:
        v2 = geo[c-1]-geo[d-1]
        r2 = np.linalg.norm(v2)
        v2 = v2/r2
        ud = np.cross(v2, cb)/(r2*m*np.linalg.norm(np.cross(v2, cb))**2)
        b_row[3*d-3] = -ud[0]
        b_row[3*d-2] = -ud[1]
        b_row[3*d-1] = -ud[2]
        vb = vb - np.dot(bc, v2)*r2*ud/(r_bc)
        vc = vc - (r_bc-r2*(np.dot(v2, cb)))*ud/r_bc

    b_row[3*b-3] = vb[0]
    b_row[3*b-2] = vb[1]
    b_row[3*b-1] = vb[2]
    b_row[3*c-3] = vc[0]
    b_row[3*c-2] = vc[1]
    b_row[3*c-1] = vc[2]

    return b_row

def wag(geo, a, b, c, d):
    """Calculates the out-of plane wag for the plane defined by atoms A, C 
    and D

    Args:
        geo (2D array): molecular geometry in Cartesian coordinates
        a (int): label of central atom A
        b (int): label of atom B
        c (int): label of atom C
        d (int): label of atom D

    Returns:
        1D array: a row of the B matrix
    """    
    b_row = np.zeros(3*len(geo))
    vab = geo[a-1]-geo[b-1]
    rab = np.linalg.norm(vab)
    eab = vab/rab
    vac = geo[a-1]-geo[c-1]
    rac = np.linalg.norm(vac)
    eac = vac/rac
    vad = geo[a-1]-geo[d-1]
    rad = np.linalg.norm(vad)
    ead = vad/rad

    normal = np.cross(eac, ead)/np.linalg.norm(np.cross(eac, ead))

    ub = normal/rab
    uc = (normal/rac) * (np.linalg.norm(np.cross(eab, ead)) /
                         np.linalg.norm(np.cross(eac, ead)))
    ud = (normal/rad) * (np.linalg.norm(np.cross(eab, eac)) /
                         np.linalg.norm(np.cross(eac, ead)))
    ua = -ub-uc-ud

    b_row[3*a-3] = ua[0]
    b_row[3*a-2] = ua[1]
    b_row[3*a-1] = ua[2]
    b_row[3*b-3] = ub[0]
    b_row[3*b-2] = ub[1]
    b_row[3*b-1] = ub[2]
    b_row[3*c-3] = uc[0]
    b_row[3*c-2] = uc[1]
    b_row[3*c-1] = uc[2]
    b_row[3*d-3] = ud[0]
    b_row[3*d-2] = ud[1]
    b_row[3*d-1] = ud[2]

    return b_row

def linear(geo, a, b, c, deg=True):
    """Calculates the angle bending internal coordinate between atoms 
    a-b-c for linear coordinates

    Args:
        geo (2D array): molecular geometry in Cartesian Coordinates
        a (int): label of atom A
        b (int): label of central atom B
        c (int): label of atom C
        deg (bool): True for degenerate states

    Returns:
        1D array: a row of the B matrix, two rows if deg = True
    """
    b_row = np.zeros(3*len(geo))
    v1 = geo[a-1]-geo[b-1]
    r1 = np.linalg.norm(v1)
    v2 = geo[c-1]-geo[b-1]
    r2 = np.linalg.norm(v2)

    u = np.array([1, 1, -(v1[0]+v1[1])/v1[2]])
    u = u/np.linalg.norm(u)     

    ua = u/r1
    uc = u/r2
    ub = -ua-uc
    b_row[3*a-3] = ua[0]
    b_row[3*a-2] = ua[1]
    b_row[3*a-1] = ua[2]
    b_row[3*b-3] = ub[0]
    b_row[3*b-2] = ub[1]
    b_row[3*b-1] = ub[2]
    b_row[3*c-3] = uc[0]
    b_row[3*c-2] = uc[1]
    b_row[3*c-1] = uc[2]
  
    if deg == True:
        b_row = np.array([b_row, np.zeros(3*len(geo))])
        v = np.cross(u, v1)
        v = v/np.linalg.norm(v)

        va = v/r1
        vc = v/r2
        vb = -va-vc
        b_row[1][3*a-3] = va[0]
        b_row[1][3*a-2] = va[1]
        b_row[1][3*a-1] = va[2]
        b_row[1][3*b-3] = vb[0]
        b_row[1][3*b-2] = vb[1]
        b_row[1][3*b-1] = vb[2]
        b_row[1][3*c-3] = vc[0]
        b_row[1][3*c-2] = vc[1]
        b_row[1][3*c-1] = vc[2]

    return b_row