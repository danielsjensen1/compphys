from itertools import repeat
from numpy import abs, argsort, empty, empty_like, exp, pi, sqrt
from scipy.integrate import simps
from scipy.linalg import eigh


class Hydrogen(object):
    def __init__(self, exponents):
        self.exponents = exponents
        self.S = empty([i for i in repeat(len(exponents), 2)])
        self.T = empty_like(self.S)
        self.A = empty_like(self.S)
        self.fill_arrays()
    
    def coulomb(self, alpha_p, alpha_q):
        return -2 * pi / (alpha_p + alpha_q)
    
    def approx_eigfunc(self, n, r):
        C = self.eigvecs[:, n]
        f = 0
        for i, alpha in enumerate(self.exponents):
            f += C[i] * exp(-alpha * r**2)
        N = 1 / sqrt(simps(abs(f)**2 * r**2, r))
        return abs(f)**2 * N**2
    
    def exact_eigfunc(self, n, r):
        return abs(2 * exp(-r))**2
    
    def fill_arrays(self):
        for i, alpha_p in enumerate(self.exponents):
            for j, alpha_q in enumerate(self.exponents):
                self.S[i, j] = self.overlap(alpha_p, alpha_q)
                self.T[i, j] = self.kinetic(alpha_p, alpha_q)
                self.A[i, j] = self.coulomb(alpha_p, alpha_q)
        self.H = self.T + self.A
    
    def kinetic(self, alpha_p, alpha_q):
        alpha_sum = alpha_p + alpha_q
        return 3 * alpha_p * alpha_q * pi**(3 / 2e0) / alpha_sum**(5 / 2e0)
    
    def overlap(self, alpha_p, alpha_q):
        return (pi / (alpha_p + alpha_q))**(3 / 2e0)
    
    def variational(self):
        eigvals, eigvecs = eigh(self.H, self.S)
        sortlist = argsort(eigvals)
        self.eigvals = eigvals[sortlist]
#        print self.eigvals
        self.eigvecs = eigvecs[sortlist]