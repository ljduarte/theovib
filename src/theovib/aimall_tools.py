from numpy import array

def get_IQA(atomicfiles, atoms):
    """Reads IQA terms from AIMALL outputs

    Args:
        atomicfiles (str): List of atomicfiles folder from AIMAll output
        atoms (array): List of atoms in the system

    Returns:
        np.array: Arrays containing E_intra, Vcl and Vxc IQA contributions
    """    
    atom_label = [atoms[i].lower() + str(i+1) for i in range(len(atoms))]
    Eintra = []
    Vcl = []
    Vxc = []
    for i in atom_label:
        f = open(atomicfiles + "/" + i + ".int")
        lines = f.readlines()
        for line in lines:
            if "E_IQA_Intra(A)" in line:
                Eintra.append(float(line.split()[2]))
    atom_int_label = []
    for i in range(len(atoms)):
        for j in range(i+1, len(atoms)):
            atom_int_label.append(
                atoms[i].lower() + str(i+1) + '_' + atoms[j].lower() + str(j+1))

    for i in atom_int_label:
        f = open(atomicfiles + "/" + i + ".int")
        lines = f.readlines()
        for line in lines:
            if "VC_IQA(A,B)" in line:
                Vcl.append(float(line.split()[5]))
            if "VX_IQA(A,B)" in line:
                Vxc.append(float(line.split()[5]))
    Eintra = array(Eintra)
    Vcl = array(Vcl)
    Vxc = array(Vxc)
    return Eintra, Vcl, Vxc

def get_energy_from_wfn(file):
    """Reads energy value from Gaussian's wavefunction file

    Args:
        file (str): file path

    Returns:
        float: energy eigenvalue from the wavefunction (Hartree)
    """
    f = open(file)
    lines = f.readlines()
    energy = float(lines[-1].split()[3])
    return energy

def get_electronic(file):
    """Reads atomic charge and dipoles from aimall output

    Args:
        file (str): .int file path

    Returns:
        charge (float): atomic charge
        dipole_x (float): x component of atomic dipole
        dipole_y (float): y component of atomic dipole
        dipole_z (float): z component of atomic dipole
    """
    f = open(file)
    lines = f.readlines()
    for line in lines:
        if "Net Charge" in line:
            charge = float(line.split()[5])
        elif "Dipole X" in line:
            dipole_x = float(line.split()[3])
        elif "Dipole Y" in line:
            dipole_y = float(line.split()[3])
        elif "Dipole Z" in line:
            dipole_z = float(line.split()[3])
    return charge, dipole_x, dipole_y, dipole_z
