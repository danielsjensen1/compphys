import matplotlib.pyplot as plt
from numpy import arange, argmin, array, linspace
from chap4.h2plus import H2Plus, Atom


H_exponents = {'STO-2G': [1.332480e0, 2.015287e-1],
               'STO-3G': [4.500225e0, 6.812745e-1, 1.513748e-1],
               'STO-4G': [1.300773e1, 1.962079e0, 4.445290e-1, 1.219492e-1],
               'STO-5G': [3.405432e1, 5.122332e0, 1.164455e0, 3.271926e-1, 1.030649e-1],
               'STO-6G': [8.294860e1, 1.245470e1, 2.833203e0, 7.999973e-1, 2.585797e-1, 8.996358e-2]}
def h2plus_calc(d, exponents):
    """
    Find the ground state of the H2+ ion.
    
    `d` - float
        Interatomic distance.
    """
    R1 = array([-d / 2e0, 0e0, 0e0])
    R2 = -R1
    Z = 1e0
    left = Atom(Z, exponents, R1)
    right = Atom(Z, exponents, R2)
    h2plus = H2Plus((left, right))
    h2plus.variational()
    return h2plus

def pes(darray, exponents):
    """
    Compute and plot the potential energy surface for H2+.
    
    `darray` - 1d array
        The (strictly increasing) array of interatomic distances. 
    """
    n = 1 #  Ground state
    energies = array([h2plus_calc(d, exponents).approx_energy(n)
                      for d in darray])
    nuc_repulsion = array([1e0 / d for d in darray])
    energies += nuc_repulsion
    print min(energies)
    print darray[argmin(energies)]
    plt.plot(darray, energies, 'r-', label='STO-4G')
    plt.xlabel(r'Interatomic distance $d$')
    plt.ylabel(r'Ground state energy (Hartrees)')
    plt.legend()
    plt.title('Potential Energy Curve for $H_{2}^{+}$')
    plt.show()
#    plt.savefig('prob01Aplot.pdf')

def plot_basis(h2plus):
    x = linspace(-3, 3, 100)
    n = 1 #  Ground state
    approx = h2plus.approx_eigfunc(n, x)
    plt.plot(x, approx, 'r-', label='STO-4G')
    plt.xlabel(r'$x$ (where $y=z=0$)')
    plt.ylabel(r'Ground state wave function')
    plt.legend()
    plt.title('Ground-state Electronic Eigenfunction of $H_{2}^{+}$ Molecule')
    plt.show()
#    plt.savefig('prob01Aplot.pdf')

if __name__ == '__main__':
    d = arange(0.5e0, 2.1e0, 0.1)
    exponents = array([13.00773, 1.962079, 0.444529, 0.1219492])
    exponents = H_exponents['STO-2G']
    pes(d, exponents)
