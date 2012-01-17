import unittest
from numpy import array, dot, exp
from numpy.linalg import inv
from chap3.prob03 import part_a


class TestPartA(unittest.TestCase):

    def test2by2(self):
        a = 2e0
        S11 = 1e0
        S12 = (1 + a) * exp(-a)
        S = array([[S11, S12], [S12, S11]])
        V = part_a(S)
        print reduce(dot, (V.transpose().conjugate(), S, V))
        print inv(V)
        print V.transpose().conjugate()


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testProb3']
    unittest.main()