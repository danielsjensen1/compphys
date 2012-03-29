from itertools import combinations, product
from numpy import dot, einsum, empty, empty_like, exp, ones, pi, sqrt
from scipy.special import erf
from chap3.prob03 import part_b as eigh
from scipy.io import savemat


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

class Helium(object):
    def __init__(self, atoms, N=1):
        """Apply the variational principle to the He atom.
        
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
        self.Q = empty((self.bsize,) * 4)
        self.fill_arrays()

    def all_orbs(self):
        """Generator that loops over all orbitals in the basis set."""
        i = 0
        for atom in self.atoms:
            for exponent in atom.exponents:
                yield i, exponent, atom.position
                i += 1

    def eigfunc(self, n, x):
        """
        Form eigenfunction `n` using LCAO.
        
        Note that currently this only does a plot of the projection of
        the eigenfunction onto the x axis (y = z = 0).
        
        `n` - positive integer
            The eigenfunction index.
        `x` - 1D array
            The x grid for the projection of the eigenfunction.
        """
        C = self.eigvecs[:, n - 1] #  Subtract 1 for 0-based indexing
        f = 0
        for i, alpha, R in self.all_orbs():
            r = x - R[0]
            f += C[i] * exp(-alpha * r ** 2)
        return f

    def energy(self, n):
        """Electronic energy only."""
        C = self.eigvecs[:, n-1] #  Subtract 1 for 0-based indexing
        E = (2 * einsum('p,q,pq', C, C, self.H) +
             einsum('ijkl,i,k,j,l', self.Q, C, C, C, C))
        return E

    def energy_tot(self, n):
        """Total energy ignoring the motion of the nuclei."""
        return self.energy(n) + self.coulomb_nuclei()

    def coulomb(self, alpha, beta, fac, R_AB_sqrd, R_P, R_C, Z_C):
        """Coulomb matrix elements."""
        R_PC = R_P - R_C
        R_PC_sqrd = dot(R_PC, R_PC)
        term1 = -2e0 * pi * Z_C / (alpha + beta)
        term2 = exp(-fac * R_AB_sqrd)
        term3 = F0((alpha + beta) * R_PC_sqrd)
        return term1 * term2 * term3

    def coulomb_nuclei(self):
        """The electrostatic energy of the nuclei."""
        combos = combinations(self.atoms, 2)
        def dist(r1, r2):
            r12 = r1 - r2
            return sqrt(dot(r12, r12))
        E = sum(atm1.charge * atm2.charge / dist(atm1.position, atm2.position)
                for atm1, atm2 in combos)
        return E

    def fill_arrays(self):
        """Fill the overlap, kinetic, and coulomb arrays."""
        for (i, alpha, R_A), (j, beta, R_B) in product(self.all_orbs(),
                                                       repeat=2):
            R_AB = R_A - R_B
            R_AB_sqrd = dot(R_AB, R_AB)
            R_P = (alpha * R_A + beta * R_B) / (alpha + beta)
            fac = alpha * beta / (alpha + beta)
            self.S[i, j] = self.overlap(alpha, beta, fac, R_AB_sqrd)
            self.T[i, j] = self.kinetic(alpha, beta, fac, R_AB_sqrd)
            self.A[i, j] = sum(self.coulomb(alpha, beta, fac, R_AB_sqrd,
                                            R_P, atom.position, atom.charge)
                               for atom in self.atoms)
        self.H = self.T + self.A
        for orbs in product(self.all_orbs(), repeat=4):
            indices, alphas, positions = zip(*orbs)
            self.Q[indices] = self.two_electron(*alphas)
#        savemat('static.mat', {'S': self.S, 'T': self.T, 'A': self.A, 'Q': self.Q})

    def kinetic(self, alpha, beta, fac, R_AB_sqrd):
        """The kinetic energy of the electrons."""
        term1 = fac * (3e0 - 2e0 * fac * R_AB_sqrd)
        term2 = (pi / (alpha + beta)) ** (1.5e0)
        term3 = exp(-fac * R_AB_sqrd)
        return term1 * term2 * term3

    def overlap(self, alpha, beta, fac, R_AB_sqrd):
        return (pi / (alpha + beta)) ** (1.5e0) * exp(-fac * R_AB_sqrd)

    def two_electron(self, alpha_p, alpha_r, alpha_q, alpha_s):
        Q = 2 * pi ** (2.5e0) / ((alpha_p + alpha_q) * (alpha_r + alpha_s) *
                                 sqrt(alpha_p + alpha_q + alpha_r + alpha_s))
        return Q

    def variational(self, Cguess=None, tol=1e-10, max_iters=100):
        if Cguess != None:
            C = Cguess #  We should check that this vector is not zero
        else:
            C = ones(self.bsize)
        #  Normalize the initial guess
        C = 1e0 / sqrt(einsum('p,pq,q', C, self.S, C)) * C
        def iteration(C):
            #  Form the Fock matrix
            F = self.H + einsum('ijkl,j,l', self.Q, C, C)
#            savemat('F.mat', {'F': F})
            self.eigvals, self.eigvecs = eigh(F, self.S)
            C = self.eigvecs[:, 0]
            E0 = self.energy(1)
            print E0, C
            return E0, C
        E0_old, C_old = iteration(C)
        diff = tol + 1e0
        iters = 1
        while (diff > tol) and (iters < max_iters):
            E0_new, C_new = iteration(C_old)
            iters += 1
            diff = abs(E0_new - E0_old)
            E0_old, C_old = E0_new, C_new
