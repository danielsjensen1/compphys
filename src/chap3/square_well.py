from numpy import argsort, array, cos, empty, pi, sin
from scipy.linalg import eigh
from sympy import diff, integrate, pprint, Symbol


class SquareWell(object):
    def __init__(self):
        self._m = Symbol('m'); m = self._m
        self._n = Symbol('n'); n = self._n
        x = Symbol('x')
        self.Smn = integrate((self.psi(m, x) * self.psi(n, x)).expand(),
                             (x, -1, 1))
        self.Hmn = integrate((-self.psi(m, x) *
                              diff(self.psi(n, x), x, 2)).expand(),
                             (x, -1, 1))
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
        return (self.Hmn).subs({self._m: m, self._n: n})
    
    def overlap(self, m, n):
        return (self.Smn).subs({self._m: m, self._n: n})
    
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
         eigvals, eigvecs = eigh(H, S)
         sortlist = argsort(eigvals)
         self.eigvals = eigvals[sortlist]
         print self.eigvals
         self.eigvecs = eigvecs[sortlist]
         