import unittest
from numpy import array, dot, exp
from numpy.linalg import inv
from scipy.linalg import eigh
from chap3.prob03 import part_a, part_b


class TestPartA(unittest.TestCase):

    def test2by2(self):
        a = 2e0
        S11 = 1e0
        S12 = (1 + a) * exp(-a)
        S = array([[S11, S12], [S12, S11]])
        V, Vdag, Vinv = part_a(S)
        print reduce(dot, (Vdag, S, V))
        print Vinv
        print inv(V)
        print Vdag
        print(dot(Vinv, V))
        print(dot(V, Vinv))

    def test5by5(self):
        print('Test 5 by 5')
        S = array([[ 1.06666667, 0., 0.15238095, 0., .05079365],
                   [ 0., 0.15238095, 0., 0.05079365, 0.],
                   [ 0.15238095, 0., 0.05079365, 0., 0.02308802],
                   [ 0., 0.05079365, 0., 0.02308802, 0.],
                   [ 0.05079365, 0., 0.02308802, 0., 0.01243201]])
        V, Vdag, Vinv = part_a(S)
        print reduce(dot, (Vdag, S, V))
        print Vinv
        print inv(V)
        print Vdag
        print(dot(Vinv, V))
        print(dot(V, Vinv))
    
class TesPartB(unittest.TestCase):
    
    def test2by2(self):
        a = 2e0
        S11 = 1e0
        S12 = (1 + a) * exp(-a)
        S = array([[S11, S12], [S12, S11]])
        H11 = -(0.5e0 + exp(-2 * a))
        H12 = -exp(-a) / 2e0 * (3 + a)
        H22 = -(0.5e0 + exp(-2 * a))
        H = array([[H11, H12], [H12, H22]])
        E, C = part_b(H, S)
        print E
        print C
        print eigh(H, S)
    


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testProb3']
    unittest.main()