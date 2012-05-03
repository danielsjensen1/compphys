from itertools import product
from numpy import (abs, append, array, complex, diff, empty, exp, gradient, 
                   invert, linspace, log, pi, r_, sin, sqrt, where)
from numpy.linalg import det, slogdet


class KronigPenny(object):
    def __init__(self, a=2e0, Delta=1e0, V0=1.5e0):
        """
        Solve the Kronig-Penny potential problem exactly or with the APW method.
        
        Parameters
        ----------
        a : float
        The spacing between barriers
        Delta : float
        The width of each barrier; must be smaller than `a`.
        V0 : float
        The height of the barriers.
        """
        self.a, self.Delta, self.V0 = (a, Delta, V0)
        #  kmax is the length of the first Brillouin zone
        self.kmax = pi / a
        if a < Delta:
            raise ValueError("`a` must be greater than `Delta`")
    
    def exact(self, E=None, Emax=30e0, numpts=2000):
        """
        Compute the exact band spectrum.
        
        Our exact expression for the band structure gives us k as a function
        of E so the user can either input a range of energies or accept the
        default range.  Note that our expression currently can only plot the
        energies that are greater than the potential `V0`.
        
        Parameters
        ----------
        E : 1D array
            Array of energies to scan.
        Emax : float
            Maximum energy to include in energy scan; not used if `E!=None`.
        numpts : int
            Number of points in energy scan; not used if `E!=None`.
        """
        a, Delta, V0 = (self.a, self.Delta, self.V0)
        if E == None:
            E = linspace(V0+1e-6, Emax, numpts)
        else:
            if all(E < V0):
                raise Exception("Invalid energy range")
        q = sqrt(2e0 * E)
        kappa = sqrt(2e0 * (E - V0))
        #  Compute the elements of the transfer matrix T
        T11 = (exp(1j * q * (a - Delta)) *
               (exp(1j * kappa * Delta) * (1e0 + kappa / q)**2 -
                exp(-1j * kappa * Delta) * (1e0 - kappa / q)**2))
        T12 = -2e0 * 1j * exp(1j * q * a) * (1e0 - kappa**2 / q**2) * sin(kappa * Delta)
        #  Solve for the eigenvalues of the 2x2 transfer matrix
        A = complex(1e0, 0e0)
        B = -2 * T11.real
        C = abs(T11)**2 - abs(T12)**2
        discriminant = B**2 - 4 * A * C
        beta = array([(-B + sign * sqrt(discriminant)) / (2 * A) 
                      for sign in (1e0, -1e0)])
        #  alpha contains both eigenvalues for each energy in array E
        alpha = q / (4 * kappa) * beta
        #  mask the unallowed energies based on the eigenvalues being complex and
        #  of magnitude 1
        mask = where( array(abs(alpha[1].real) <= 1e0) * array(discriminant <= 0e0))
        k = log(alpha) / (1j * a)
        #  We know that the final answer must be all real so we drop the
        #  imaginary part.  We could do a check here (using the mask) to assert
        #  that the imaginary part really is negligible.
        k = k.real
        kexact, Eexact = (append(k[0][mask], k[1][mask]), append(E[mask], E[mask]))
        return kexact, Eexact
    
    def fill_arrays(self, k, E):
        """
        Fill the Hamiltonian (H) and overlap (S) matrices.
        
        Parameters
        ----------
        k : float
            The k point of interest.
        E : float
            The trial energy.
        """
        a, V0, Delta = (self.a, self.V0, self.Delta)
        kappa = sqrt(2 * (E - V0))
        def q(i):
            return k + 2 * pi * i / a
        def C(i):
            return sin((kappa + q(i)) * Delta / 2e0) / sin(kappa * Delta)
        def D(i):
            return sin((kappa - q(i)) * Delta / 2e0) / sin(kappa * Delta)
        S = empty((self.bsize,) * 2)
        H = empty((self.bsize,) * 2)
        for m, n in product(self.m, repeat=2):
            qm, qn = (q(m), q(n))
            qmn = qm - qn
            Cm, Cn = (C(m), C(n))
            Dm, Dn = (D(m), D(n))
            if m == n:
                Sext = a - Delta
            else:
                Sext = -2e0 / qmn * sin(qmn * Delta / 2e0)
            term1 = (Cm * Cn + Dm * Dn) * Delta
            term2 = (Cm * Dn + Dm * Cn) * sin(kappa * Delta) / kappa
            Sint = term1 + term2
            S[m, n] = Sext + Sint
            Hext = 0.5e0 * qm * qn * Sext
            Hint = kappa**2 / 2e0 * (term1 - term2) + V0 * Sint
#            Hint = E * Sint
            H[m, n] = Hext + Hint
        return H, S
    
    def find_roots(self, Earray, sign, logdet):
        """
        Find the roots of the natural logarithm of some quantity.
        
        Parameters
        ----------
        Earray : 1D array
            Energy scan array.
        sign : 1D array
            The sign of the determinant for each energy.
        logdet : 1D array
            The natural logarithm of the absolute value of the determinant for
            each energy.
        """
        #  Find the indices immediately to the left of the crossings
        crossings = where(diff(sign))[0]
        #  The exact zero is somewhere between the left and right indices
        energies = (Earray[crossings] + Earray[crossings+1]) / 2e0
        slope = gradient(logdet)
        mask = where(invert(slope[crossings-1] > 0e0))
        return energies[mask]       
        
    def find_roots_old(self, Earray, sign, logdet):
        """
        A more primitive root finder that simply locates the local minima.
        """
        #  Find all local minima
        indices = r_[False, logdet[1:] < logdet[:-1]] & \
                  r_[logdet[:-1] < logdet[1:], False] & \
                  r_[False, diff(sign) != 0]
        #  The first index has to be false since it is being compared with the
        #  last index
        return Earray[indices]

    def apw(self, karray, Earray, m=range(-6,7)):
        x = []
        y = []
        for k in karray:
            sign, logdet = self.scan_energies(k, Earray, m=m)
            energies = self.find_roots(Earray, sign, logdet)
            map(x.append, [k]*len(energies))
            map(y.append, energies)
        return x, y

    def scan_energies(self, k, Earray, m=range(-6,7)):
        self.bsize = len(m)
        self.m = m
        sign = []
        logdet = []
        for E in Earray:
            H, S = self.fill_arrays(k, E)
#            val = det(H - E * S)
            val = slogdet(H-E*S)
            sign.append(val[0])
            logdet.append(val[1])
        return array(sign), array(logdet)
