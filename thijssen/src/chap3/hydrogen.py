from numpy import empty, empty_like, exp, pi, sqrt
from chap3.prob03 import part_b as eigh


class Hydrogen(object):
    def __init__(self, exponents):
        self.exponents = exponents
        self.S = empty((len(exponents),) * 2)
        self.T = empty_like(self.S)
        self.A = empty_like(self.S)
        self.fill_arrays()
    
    def coulomb(self, alpha_p, alpha_q):
        return -2 * pi / (alpha_p + alpha_q)
    
    def approx_eigfunc(self, n, r):
        C = self.eigvecs[:, n - 1] #  Subtract 1 for 0-based indexing
        f = 0
        #  Form the eigenfunction using the correct linear combination 
        #  of the basis functions.
        for i, alpha in enumerate(self.exponents):
            f += C[i] * exp(-alpha * r**2)
        return f
    
    def approx_energy(self, n):
        return self.eigvals[n - 1] #  Subtract 1 for 0-based indexing 
    
    def exact_eigfunc(self, n, r):
        #  Currently only the ground state eigenfunction is implemented
        return exp(-r) / sqrt(pi)
    
    def exact_energy(self, n):
        return -1e0 / float(2 * n**2)
    
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
        self.eigvals, self.eigvecs = eigh(self.H, self.S)
