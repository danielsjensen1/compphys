from itertools import product
from numpy import (abs, append, array, complex, diff, empty, exp, gradient, linspace, log, pi, 
                   sign, sin, sqrt, where)
from scipy.linalg import det
from numpy.linalg import slogdet
import matplotlib.pyplot as plt


def plot_exact(a=2e0, Delta=1e0, V0=1.5e0, Vmax=30e0, numpts=500):
    E = linspace(V0+1e-6, Vmax, numpts)
    q = sqrt(2e0 * E)
    kappa = sqrt(2e0 * (E - V0))
    T11 = (exp(1j * q * (a - Delta)) *
           (exp(1j * kappa * Delta) * (1e0 + kappa / q)**2) -
           (exp(-1j * kappa * Delta) * (1e0 - kappa / q)**2))
    T12 = -2e0 * 1j * exp(1j * q * a) * (1e0 - kappa**2 / q**2) * sin(kappa * Delta)
#    determinant = (abs(T11)**2-abs(T12)**2)*(q/(4*kappa))**2
#    print(determinant)
#    mask = where(abs(determinant - 1e0)<1e-1)
#    print('T11=', T11)
#    print('T12=', T12)
    A = complex(1e0, 0e0)
#    b = -2 * T11.real
    B = -T11 - T11.conjugate()
    C = abs(T11)**2 - abs(T12)**2
    beta = array([(-B + sign * sqrt(B**2 - 4 * A * C)) / (2 * A) 
                  for sign in (1e0, -1e0)])
    alpha = q / (4 * kappa) * beta
    mask = where(abs(alpha[0].real) <= 1e0)
#    mask = where(abs(sum(alpha)/2e0) <= 1e0+1e-6)
    print(abs(alpha))
#    print('alpha =', alpha)
    k = log(alpha) / (1j * a)
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    kexact, Eexact = (append(k[0][mask], k[1][mask]), append(E[mask], E[mask]))
    ax1.plot(kexact, Eexact, 'b.')
    ax1.set_title('Exact Band Structure')
    ax1.set_xlabel('$k$ point')
    ax1.set_ylabel('Energy')
#    kmax = 2 * pi / a
#    E = linspace(0, Vmax, numpts)
#    kplus = (sqrt(2 * E) - kmax / 2e0) % kmax - kmax/2e0
#    kminus = -kplus
#    ax1.plot(kplus, E, 'r-', kminus, E, 'r-')
#    ax1.set_xlim((0, kmax/2e0))
    plt.show()
    return kexact, Eexact

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
#            print(E, slogdet(H-E*S))
#            val = det(H - E * S)
            val = slogdet(H-E*S)
            vals.append(val)
#            print('{0:.5e}, {1:.5g}'.format(E, val))
        vals = array(vals)
        sign, logdet = vals[:, 0], vals[:, 1]
        slope = gradient(logdet)
        crossings = where(diff(sign))
        energies = Earray[crossings]
        print slope[crossings]
        allowed_energies = energies[where(slope[crossings[0]-1] < 0e0)]
#        allowed_energies = energies[where(slope[crossings-1] < 0e0 and
#                                          slope[crossings+2] > 0e0)]
        print(allowed_energies)
#        print slope[crossings]
#        for crossing in crossings[0]:
#            print(slope[crossing:crossing+2], Earray[crossing])
#        print('crossings=', Earray[crossings])
        plt.plot(Earray, sign, label='sign')
        plt.plot(Earray, logdet, label='log(det)')
        plt.title("log(Determinant) vs. Energy")
        plt.xlabel('energy')
        plt.ylabel('logdet')
        plt.show()
        return allowed_energies


if __name__ == '__main__':
    a, Delta, V0, Vmax = (2e0, 1e0, 1.5e0, 30e0)
    kexact, Eexact = plot_exact(a, Delta, V0, Vmax)
    apw = APW(a, Delta, V0)
    Earray = linspace(apw.V0+1e-6, 30e0, 200)
    k = 0.1e0
    apw.scan_energies(k, Earray)
#    plt.plot(kexact, Eexact, 'b.')
#    for k in linspace(-pi/2e0, pi/2e0, 30):
#        energies = apw.scan_energies(k, Earray)
#        plt.plot([k]*len(energies), energies, 'r.')
#    plt.title('Exact and APW Band Structure')
#    plt.xlabel('$k$ point')
#    plt.ylabel('Energy')
#    plt.show()