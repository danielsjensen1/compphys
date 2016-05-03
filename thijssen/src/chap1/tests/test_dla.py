import unittest
from chap1.dla import DLA


class TestDLA(unittest.TestCase):


    def test5x5(self):
        N = 11
        dla = DLA(size=(N, N), occupied=None)
        print(dla.grid)
        dla.walker((0, 0))
        print(dla.grid)
        dla.cluster((0, 0), Nsites=5)
        print(dla.grid)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()