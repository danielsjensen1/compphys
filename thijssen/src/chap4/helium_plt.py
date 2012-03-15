import matplotlib.pyplot as plt
from numpy import array, linspace
from chap4.helium import Helium, Atom


#  Place the atoms one bohr radius apart
R1 = array([0e0, 0e0, 0e0])
Z = 2e0
exponents = array([0.298073e0, 1.242567e0, 5.782948e0, 38.474970e0])
he1 = Atom(Z, exponents, R1)
helium = Helium((he1,))
helium.variational()
def print_energies():
    for n in range(1, helium.bsize + 1):
        approx = helium.energy(n)
        print('E{0} = {1:+.7e}'.format(n, approx))
def plot_groundstate():
    x = linspace(-3, 3, 100)
    n = 1 #  Ground state
    approx = helium.eigfunc(n, x)
    plt.plot(x, approx, 'x-', label='STO-4G')
    plt.xlabel(r'$x$ (where $y=z=0$)')
    plt.ylabel(r'Ground state wave function')
    plt.legend()
    plt.title('Ground-state Electronic Eigenfunction of $H_{2}^{+}$ Molecule')
    plt.show()
#    plt.savefig('prob01Aplot.pdf')

if __name__ == '__main__':
    print_energies()
#    plot_groundstate()
