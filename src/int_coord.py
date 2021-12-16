import numpy as np 

def bond(geo, a, b):
    """calculate the bond internal coordinate between a and b in the molecule

    :param geo: molecule geometry in Cartesian Coordinate
    :type geo: 2D array  
    :param a: label of atom a
    :type a: int
    :param b: label of atom b
    :type b: int]
    :return: b_row
    :rtype: 1D array 
    """
    v = geo[a] - geo[b]
    v = v/np.linalg.norm(v)
    b_row = np.zeros(3*len(geo))
    b_row[3*a-3] = v[0]
    b_row[3*a-2] = v[1]
    b_row[3*a-1] = v[2]
    b_row[3*b-3] = -v[0]
    b_row[3*b-2] = -v[1]
    b_row[3*b-1] = -v[2]

    return b_row

def angle(geo, a, b, c):
    """[summary]

    :param geo: [description]
    :type geo: [type]
    :param a: [description]
    :type a: [type]
    :param b: [description]
    :type b: [type]
    :param c: [description]
    :type c: [type]
    :return: [description]
    :rtype: [type]
    """

    v1 = geo[a]-geo[b]
    v2 = geo[c]-geo[b]
    normal = np.cross(v1,v2)
    u1 = np.cross(normal, v1)
    u1 = u1/np.linalg.norm(u1)
    u2 = np.cross(normal, v2)
    u2 = u1/np.linalg.norm(u2)

    b_row = np.zeros(3*len(geo))
    b_row[3*a-3] = u1[0]
    b_row[3*a-2] = u1[1]
    b_row[3*a-1] = u1[2]
    b_row[3*b-3] = u2[0]
    b_row[3*b-2] = u2[1]
    b_row[3*b-1] = u2[2]
    b_row[3*c-3] = -(u1[0]+u2[0])
    b_row[3*c-2] = -(u1[1]+u2[1])
    b_row[3*c-1] = -(u1[2]+u2[2])

    return b_row

def dihedral(geo, a, b, c, d):
    pass

