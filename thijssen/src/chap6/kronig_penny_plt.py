from math import pi
import time
import matplotlib.pyplot as plt
from numpy import linspace
from chap6.kronig_penny import KronigPenny


a, Delta, V0 = (2e0, 1.5e0, 2e0)
#Emax = 40 * V0
Emax = 20e0
kp = KronigPenny(a, Delta, V0)

def plot_apw(kpts=30, Epts=400, m=range(-4,5)):
    Earray = linspace(V0+1e-6, Emax, Epts)
    karray = linspace(-kp.kmax, kp.kmax, kpts)
    tic = time.clock()
    x, y = kp.apw(karray, Earray, m=m)
    toc = time.clock()
    print('time = {0} for m = {1}'.format(toc - tic, len(m)))
    plt.plot(x, y, 'ro', label='APW')
    kexact, Eexact = kp.exact(numpts=20*kpts, Emax=Emax)
    plt.plot(kexact, Eexact, 'b.', ms=2.5, label='Exact')
    plt.title('Exact and APW Band Structure ({0} plane waves)'.format(len(m)))
    plt.xlabel('$k$ point')
    plt.ylabel('Energy')
    plt.legend()
    plt.show()
#    name = 'APW_m{0}_a2_D1_V2.pdf'.format(len(m))
#    plt.savefig(name)
#    plt.clf()

def plot_exact(numpts=2000):
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    kexact, Eexact = kp.exact(numpts=numpts, Emax=Emax)
    ax1.plot(kexact, Eexact, 'b.', ms=2.5)
    ax1.hlines(y=V0, xmin=-kp.kmax, xmax=kp.kmax, label='$V_0$',
                linewidth=2.5, color='grey', linestyle='--')
    ax1.set_ylim((0, Emax))
    ax1.set_title('Exact Band Structure ($\Delta={0:.3f}$)'.format(kp.Delta))
    ax1.set_xlabel('$k$ point')
    ax1.set_ylabel('Energy')
    ax1.legend()
#    name = 'Exact_a2_D{0:.3f}_V2.pdf'.format(kp.Delta)
#    plt.savefig(name)
#    plt.clf()
    plt.show()

def plot_logdet(exact=True):
    Earray = linspace(V0+1e-6, Emax, 200)
    k = 0.1625
    sign, logdet = kp.scan_energies(k, Earray, m=range(-3,4))
    determinant = r"\left(\left|H-SE\right|\right)"
    plt.plot(Earray, sign, 'r-', label='sign${0}$'.format(determinant))
    plt.plot(Earray, logdet, 'b-', lw=2, label=r'$\log{0}$'.format(determinant))
    plt.title(r"$\log\left(\left|H-SE\right|\right)$ vs. $E$")
    plt.xlabel('energy')
    plt.ylabel('logdet')
    plt.legend()
    energies = kp.find_roots(Earray, sign, logdet)
    print(energies)
    plt.show()

if __name__ == '__main__':
    plot_logdet()
    for Delta in linspace(0, a, 8)[1:-1]:
        kp.Delta = Delta
        plot_exact()
    plot_apw(m=range(-3,3+1))
#    for M in range(1,7):
#        plot_apw(m=range(-M,M+1))