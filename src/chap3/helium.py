from itertools import product, repeat
from numpy import argsort, empty, empty_like, ndindex, ones, pi, sqrt
from scipy.linalg import eigh
from chap3.hydrogen import Hydrogen


class Helium(Hydrogen):
    def __init__(self, exponents):
        self.Q = empty([i for i in repeat(len(exponents), 4)])
        Hydrogen.__init__(self, exponents)
        
    
    def coulomb(self, alpha_p, alpha_q):
        return 2 * Hydrogen.coulomb(self, alpha_p, alpha_q)
    
    def fill_arrays(self):
        for index in ndindex(self.Q.shape):
            alphas = [self.exponents[i] for i in index]
            self.Q[index] = self.hartree(*alphas)
        Hydrogen.fill_arrays(self)
    
    def fock_and_solve(self, C):
        F = empty([i for i in repeat(len(exponents), 2)])
        Qpq = empty_like(F)
        for index in ndindex(Qpq.shape):
            p, q = index
            Qpq[index] = sum([self.Q[p, r, q, s] * C[r] * C[s] for (r, s)
                              in product(range(len(self.exponents)), repeat=2)])
        self.F = self.H + Qpq
        eigvals, eigvecs = eigh(self.F, self.S)
        #  eigh returns normalized eigenvectors so we don't normalize them again
        sortlist = argsort(eigvals)
        C = eigvecs[sortlist][:, 0]
        Eg = sum([2 * C[p] * C[q] * self.H[p, q] for (p, q) in
                  product(range(len(self.exponents)), repeat=2)])
        Eg += sum([self.Q[p, r, q, s] * C[p] * C[q] * C[r] * C[s] for 
                   (p, q, r, s) in product(range(len(self.exponents)), repeat=4)])
#        self.eigvals = eigvals[sortlist]
#        self.eigvecs = eigvecs[sortlist]
#        print self.eigvals
        return Eg, C
#        return eigvals[sortlist][0], eigvecs[sortlist][:, 0]
    
    def hartree(self, alpha_p, alpha_r, alpha_q, alpha_s):
        """Notice that the order is prqs as in Thijssen."""
        num = 2 * pi**(5 / 2e0)
        den = ((alpha_p + alpha_q) * (alpha_r + alpha_s) *
               sqrt(alpha_p + alpha_q + alpha_r + alpha_s))
        return num / den
    
#    def kinetic(self, alpha_p, alpha_q):
#        return 2 * Hydrogen.coulomb(self, alpha_p, alpha_q)
    
    def self_consistency(self, C=None, maxiters=100, tol=1e-14):
        if C is None:
            C = ones(len(self.exponents))
        N = 1e0 / sqrt(sum([C[p] * self.S[p, q] * C[q] for (p,q) in
                          product(range(len(self.exponents)), repeat=2)]))
        C = N * C
        
        #  Initialize Eold and Cold
        Eold, Cold = self.fock_and_solve(C)
        for i in range(maxiters):
            Enew, Cnew = self.fock_and_solve(Cold)
            print Enew
            if abs(Enew - Eold) < tol:
                break
            Eold, Cold = Enew, Cnew
#        C = Cnew
#        N = 1e0 / sqrt(sum([C[p] * self.S[p, q] * C[q] for (p,q) in
#                          product(range(len(self.exponents)), repeat=2)]))
#        C = N * C
        
#        print Eg
            
        

if __name__ == '__main__':
    exponents = [0.298073, 1.242567, 5.782948, 38.474970]
    helium = Helium(exponents)
    helium.self_consistency()