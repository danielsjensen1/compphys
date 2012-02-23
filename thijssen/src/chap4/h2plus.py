from itertools import product
from numpy import dot, empty, empty_like, exp, pi, sqrt
from scipy.special import erf
from chap3.prob03 import part_b as eigh


class Atom(object):
    def __init__(self, charge, exponents, position):
        """
        Store the basis information for each atom
        
        `charge` - integer
            Nuclear charge of the atom.
        `exponents` - 1D array of floats
            The exponents determining the width of the Gaussian basis functions.
        `position` - 1x3 array of floats
            The location of the atom in Cartesian coordinates.
        """
        self.charge = charge
        self.exponents = exponents
        self.position = position
        
def F0(t):
    """F0 function described on page 75 of Thijssen equation 4.114."""
    if t < 1e-6: #  limit as t->0 of F0 is 1e0
        return 1e0
    else:
        thalf = sqrt(t)
        return sqrt(pi) / (2e0 * thalf) * erf(sqrt(t))

class H2Plus(object):
    def __init__(self, atoms, N=1):
        """Apply the variational principle to the H2+ molecule.
        
        `atoms` - list or tuple of class `Atom`
            All of the atoms in the molecule and their properties.
        `N` - integer
            Total number of electrons.
        """
        self.atoms = atoms
        self.N = N
        self.bsize = sum(len(atom.exponents) for atom in atoms)
        self.S = empty((self.bsize,) * 2)
        self.T = empty_like(self.S)
        self.A = empty_like(self.S)
        self.fill_arrays()
    
    def all_orbs(self):
        """Generator that loops over all orbitals in the basis set."""
        i = 0
        for atom in self.atoms:
            for exponent in atom.exponents:
                yield i, exponent, atom.position
                i += 1
    
    def coulomb(self, alpha, beta, fac, R_AB_sqrd, R_P, R_C, Z_C):
        """Coulomb matrix elements."""
        R_PC = R_P - R_C
        R_PC_sqrd = dot(R_PC, R_PC)
        term1 = -2e0 * pi * Z_C / (alpha + beta)
        term2 = exp(-fac * R_AB_sqrd)
        term3 = F0((alpha + beta) * R_PC_sqrd)
        return term1 * term2 * term3
    
    def approx_eigfunc(self, n, x):
        """
        Form eigenfunction `n` using LCAO.
        
        Note that currently this only does a plot of the projection of
        the eigenfunction with y=z=0.
        
        `n` - positive integer
            The eigenfunction index.
        `x` - 1D array
            The x grid for the projection of the eigenfunction.
        """
        C = self.eigvecs[:, n - 1] #  Subtract 1 for 0-based indexing
        f = 0
        for i, alpha, R in self.all_orbs():
            r = x - R[0]
            f += C[i] * exp(-alpha * r**2)
        return f
    
    def approx_energy(self, n):
        return self.eigvals[n - 1] #  Subtract 1 for 0-based indexing 
    
    def fill_arrays(self):
        for (i, alpha, R_A), (j, beta, R_B) in product(self.all_orbs(), 
                                                       repeat=2):
            R_AB = R_A - R_B
            R_AB_sqrd = dot(R_AB, R_AB)
            R_P = (alpha  * R_A + beta * R_B) / (alpha + beta)
            fac = alpha * beta / (alpha + beta)
            self.S[i, j] = self.overlap(alpha, beta, fac, R_AB_sqrd)
            self.T[i, j] = self.kinetic(alpha, beta, fac, R_AB_sqrd)
            self.A[i, j] = sum(self.coulomb(alpha, beta, fac, R_AB_sqrd, 
                                            R_P, atom.position, atom.charge)
                               for atom in self.atoms)
        self.H = self.T + self.A
    
    def kinetic(self, alpha, beta, fac, R_AB_sqrd):
        term1 = fac * (3e0 - 2e0 * fac * R_AB_sqrd)
        term2 = (pi / (alpha + beta))**(1.5e0)
        term3 = exp(-fac * R_AB_sqrd)
        return term1 * term2 * term3
    
    def overlap(self, alpha, beta, fac, R_AB_sqrd):
        return (pi / (alpha + beta))**(1.5e0) * exp(-fac * R_AB_sqrd) 
    
    def variational(self):
        self.eigvals, self.eigvecs = eigh(self.H, self.S)
