import unittest

from numpy import linalg
from internal import *
from molecule import Molecule

water = Molecule(['O', 'H', 'H'], np.array([[0.000000,  0.000000,  0.004316],
                                            [0.000000, -0.763369, -0.580667],
                                            [0.000000,  0.763369, -0.580667]]))

ethylene = Molecule(['C', 'C', 'H', 'H', 'H', 'H'],
                    np.array([[-0.666350,  0.000000, 0.000000],
                              [0.666350,  0.000000, 0.0000000],
                              [-1.244722,  0.916769, 0.000000],
                              [-1.244722, -0.916769, 0.000000],
                              [1.244722,  0.916769, 0.000000],
                              [1.244722, -0.916769, 0.000000]]))

co2 = Molecule(['C', 'O', 'O'], np.array([[0.000000, 0.000000,  0.000000],
                                          [0.000000, 0.000000, -1.160431],
                                          [0.000000, 0.000000,  1.16043]]))


class TestMoleculeClass(unittest.TestCase):
    def test_read_gaussian(self):
        molecule = Molecule.read_gaussian('./code/molecule.gjf')
        for pair in zip(water.positions[0], molecule.positions[0]):
            self.assertAlmostEqual(pair[0], pair[1], 4)


class TestInternalCoords(unittest.TestCase):

    def test_bond(self):
        result = bond(water.positions, 2, 1)
        expected_list = [0.0, 0.7937401916095669, 0.608256942602062,
                         0.0, -0.7937401916095669, -0.608256942602062, 0.0, 0.0, 0.0]
        for pair in zip(expected_list, result):
            self.assertAlmostEqual(pair[0], pair[1], 4)

    def test_angle(self):
        result = angle(water.positions, 2, 1, 3)
        expected_list = [0.0, 0.0, -1.6506336257399439, 0.0, -0.6324546578098675,
                         0.8253168128699719, 0.0, 0.6324546578098675, 0.8253168128699719]
        for pair in zip(expected_list, result):
            self.assertAlmostEqual(pair[0], pair[1], 4)

    def test_torsion(self):
        result = torsion(ethylene.positions, [3, 4], 1, 2, [5, 6])
        expected_list = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.545394,
                         0.0, 0.0, 0.545394, 0.0, 0.0, 0.545394, 0.0, 0.0, -0.545394]
        for pair in zip(expected_list, result):
            self.assertAlmostEqual(pair[0], pair[1], 4)

    def test_wag(self):
        result = wag(ethylene.positions, 2, 1, 6, 5)
        expected_list = [0.0, 0.0, 0.750356, 0.0, 0.0, -2.479150, 0.0,
                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.864397, 0.0, 0.0, 0.864397]
        for pair in zip(expected_list, result):
            self.assertAlmostEqual(pair[0], pair[1], 3)

    def test_linear(self):
        vectors = linear(co2.positions, 2, 1, 3, deg=True)
        result = np.dot(vectors[0][0:3], vectors[1][0:3])
        self.assertAlmostEqual(result, 0, 5)


if __name__ == '__main__':
    unittest.main()
