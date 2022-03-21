
def find(str, ch):
    indexes = []
    line = str.split()
    for i in range(len(line)):
        if line[i] == ch:
            indexes.append(i)
    return indexes


class Input:
    def __init__(self):
        self.molecule = None
        self.folder = None
        self.bond = []
        self.angle = []
        self.linear_angle = []
        self.torsion = []
        self.oop_wag = []
        self.delta = None

    @classmethod
    def read_text(cls, file):
        a = cls()
        f = open(file, 'r').readlines()
        for i in range(len(f)):
            if f[i] == 'MOLECULE:\n':
                a.molecule = f[i+1]
            elif f[i] == 'FOLDER:\n':
                a.folder = f[i+1].rstrip('\n')
            elif f[i] == 'BOND:\n':
                bond = []
                j = i+1
                while f[j] != "---\n":
                    bond.append([int(f[j].split()[0]),
                                 int(f[j].split()[1])])
                    j = j+1
                a.bond = bond
            elif f[i] == 'ANGLE:\n':
                angle = []
                j = i+1
                while f[j] != "---\n":
                    angle.append([int(f[j].split()[0]), int(
                        f[j].split()[1]), int(f[j].split()[2])])
                    j = j+1
                a.angle = angle
            elif f[i] == 'LINEAR ANGLE:\n':
                linear_angle = []
                j = i+1
                while f[j] != "---\n":
                    linear_angle.append([int(f[j].split()[0]), int(
                        f[j].split()[1]), int(f[j].split()[2])])
                    j = j+1
                a.linear_angle = linear_angle
            elif f[i] == 'TORSION:\n':
                torsion = []
                j = i+1
                while f[j] != "---\n":
                    sep = find(f[j], '-')
                    print(f[j])
                    n = f[j].split()[:sep[0]]
                    n = list(map(int, n))
                    m = f[j].split()[sep[1]+1:]
                    m = list(map(int, m))
                    torsion.append(
                        [n, int(f[j].split()[sep[0]+1]), int(f[j].split()[sep[1]-1]), m])
                    j = j+1
                a.torsion = torsion
                print(a.torsion)
            elif f[i] == 'OOP WAG:\n':
                oop_wag = []
                j = i+1
                while f[j] != "---\n":
                    oop_wag.append([int(f[j].split()[0]), int(
                        f[j].split()[1]), int(f[j].split()[2]), int(f[j].split()[3])])
                    j = j+1
                a.oop_wag = oop_wag
            elif f[i] == "DELTA:\n":
                a.delta = float(f[i+1])
        return a
