from numpy import cos, empty, pi, sin
#from scipy.linalg import eigh
from chap3.prob03 import part_b as eigh
from sympy import diff, integrate, Symbol


class SquareWell(object):
    def __init__(self):
        self._m = Symbol('m'); m = self._m
        self._n = Symbol('n'); n = self._n
        x = Symbol('x')
        def Smn(m, n):
            return integrate((self.psi(m, x) * self.psi(n, x)).expand(basic=True),
                             (x, -1, 1))
        def Hmn(m, n):
            return integrate((-self.psi(m, x) *
                              diff(self.psi(n, x), x, 2)).expand(basic=True),
                             (x, -1, 1))
        self.Hmn, self.Smn = Hmn, Smn
    def getwfs(self):
        'Number of wave functions used in approximation.'
        return self.__wfs
    def setwfs(self, value):
        self.__wfs = value
    wfs = property(getwfs, setwfs)
    
    def approx_eigfunc(self, n, x):
        C = self.eigvecs[:, n]
        f = 0
        for i in range(0, self.wfs):
            f += C[i] * self.psi(i, x)
        return f
    
    def exact_eigfunc(self, n, x):
        kn = n * pi / 2e0
        if n % 2 == 1:
            return cos(kn * x)
        else:
            return sin(kn * x)
    
    def kinetic(self, m, n):
        return self.Hmn(m, n)
    
    def overlap(self, m, n):
        return self.Smn(m, n)
    
    def psi(self, n, x):
        return x**n * (x - 1) * (x + 1)
    
    def variational(self, wfs=5):
        self.wfs = wfs
        H = empty([wfs, wfs])
        S = empty([wfs, wfs])
        for m in range(0, wfs):
            for n in range(0, wfs):
                H[m, n] = (self.kinetic(m, n)).evalf(15)
                S[m, n] = (self.overlap(m, n)).evalf(15)
#         print('S')
#         print(S)
#         print('H')
#         print(H)
        self.eigvals, self.eigvecs = eigh(H, S)
#        sortlist = argsort(eigvals)
#        self.eigvals = eigvals[sortlist]
        print self.eigvals
        print self.eigvecs
#        self.eigvecs = eigvecs[sortlist]