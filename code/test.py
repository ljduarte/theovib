import unittest
from int_coord import *
from molecule import Molecule

water = Molecule(['O', 'H', 'H'], np.array([[0.000000,  0.000000,  0.004316],
                                            [0.000000, -0.763369, -0.580667],
                                            [0.000000,  0.763369, -0.580667]]))

ethylene = Molecule(['C', 'C', 'H', 'H', 'H', 'H'], np.array([[-0.666350,  0.000000, 0.000000],
                                                              [ 0.666350,  0.000000, 0.0000000],
                                                              [-1.244722,  0.916769, 0.000000],
                                                              [-1.244722, -0.916769, 0.000000],
                                                              [ 1.244722,  0.916769, 0.000000],
                                                              [ 1.244722, -0.916769, 0.000000]]))

class TestInternalCoords(unittest.TestCase):

    def test_bond(self):
        result = bond(water.positions, 2, 1)
        expected_list = [0.0, 0.7937401916095669, 0.608256942602062,
                         0.0, -0.7937401916095669, -0.608256942602062, 0.0, 0.0, 0.0]
        for pair in zip(expected_list, result):
            self.assertAlmostEqual(pair[0], pair[1], 4)

    def test_angle(self):
        result = angle(water.positions, 2, 1, 3)
        print(result)
        expected_list = [0.0, 0.0, -1.6506336257399439, 0.0, -0.6324546578098675,
                         0.8253168128699719, 0.0, 0.6324546578098675, 0.8253168128699719]
        for pair in zip(expected_list, result):
            self.assertAlmostEqual(pair[0], pair[1], 4)
    
    def test_torsion(self):
        result = torsion(ethylene.positions, [3, 4], 1, 2, [5, 6])
        print(result)
        expected_list = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.545394, 0.0, 0.0, 0.545394, 0.0, 0.0, 0.545394, 0.0, 0.0, -0.545394]
        for pair in zip(expected_list, result):
            self.assertAlmostEqual(pair[0], pair[1], 4)


if __name__ == '__main__':
    unittest.main()
