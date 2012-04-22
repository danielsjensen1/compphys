from itertools import product
from numpy import (abs, array, complex, diff, empty, exp, linspace, log, pi, 
                   piecewise, sign, sin, sqrt, where)
from scipy.linalg import det
from numpy.linalg import slogdet
import matplotlib.pyplot as plt


def plot_potential(a=2e0, Delta=1e0, V0=1.5e0, barriers=4, numpts=1000):
    """
    Plot the Kronig-Penny potential
    
    Parameters
    ----------
    a : real
        The direct-space lattice vector
    Delta : real
        The width of the potential.  This width must be less than `a`.
    V0 : real
        The height of the barrier
    """
    x = linspace(0e0, barriers*a, num=numpts)
    y = piecewise(x, [(x - Delta / 2e0) % a < Delta], [V0])
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    ax1.plot(x, y)
    ax1.set_title('Kronig-Penny Potential')
    ax1.set_xticks(linspace(a/2e0, a/2e0+(barriers-1)*a, barriers))
    ax1.set_xticklabels(['${0}a$'.format(i) for i in range(barriers)])
    ax1.set_yticks([V0])
    ax1.set_yticklabels(['$V_0$'])
    ax1.set_ylim((0, V0+V0/2e0))
    plt.show()

def plot_exact(a=2e0, Delta=1e0, V0=1.5e0, Vmax=30e0):
    E = linspace(V0+1e-6, Vmax, 1000)
    q = sqrt(2e0 * E)
    kappa = sqrt(2e0 * (E - V0))
    T11 = (exp(1e0j * q * (a - Delta)) *
           (exp(1e0j * kappa * Delta) * (1e0 + kappa / q)**2) -
           (exp(-1e0j * kappa * Delta) * (1e0 - kappa / q)**2))
    T12 = -2e0 * 1e0j * exp(1j * q * a) * (1e0 - kappa**2 / q**2)
    determinant = (abs(T11)**2-abs(T12)**2)*(q/(4*kappa))**2
    print(determinant)
    mask = where(abs(determinant - 1e0)<1e-2)
#    print('T11=', T11)
#    print('T12=', T12)
    A = complex(1e0, 0e0)
#    b = -2 * T11.real
    B = -T11 - T11.conjugate()
    C = abs(T11)**2 - abs(T12)**2
    beta = array([(-B + sign * sqrt(B**2 - 4 * A * C)) / (2 * A) 
                  for sign in (1e0, -1e0)])
    alpha = q / (4 * kappa) * beta
    print(abs(alpha))
#    print('alpha =', alpha)
    k = log(alpha) / (1j * a)
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
#    ax1.plot(k[0], E, k[1], E)
    ax1.plot(k[0][mask], E[mask], k[1][mask], E[mask])
    k = linspace(-pi/a, pi/a, len(E[mask]))
    E = k**2 / 2e0
    ax1.plot(k, E)
    print(k[0])
    print(pi/a)
    plt.show()

class APW(object):
    def __init__(self, a=2e0, Delta=1e0, V0=1.5e0, m=range(-6,7)):
        #  Should I exclude m = 0?
        self.a, self.Delta, self.V0, self.m = (a, Delta, V0, m)
        if a < Delta:
            raise ValueError('`a` must be greater than `Delta`')
        self.bsize = len(m)
    
    def fill_arrays(self, k, E):
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
    
    def scan_energies(self, k, Earray):
        vals = []
        for E in Earray:
            H, S = self.fill_arrays(k, E)
            print(E, slogdet(H-E*S))
            val = det(H - E * S)
            vals.append(val)
#            print('{0:.5e}, {1:.5g}'.format(E, val))
        vals = array(vals)
        print('crossings=', Earray[where(diff(sign(vals)))])
        plt.plot(Earray, vals)
        plt.show()


if __name__ == '__main__':
    apw = APW()
    k = 0.9e0
    Earray = linspace(apw.V0+1e-6, 30e0, 100)
    apw.scan_energies(k, Earray)
#    plot_potential()
    plot_exact()