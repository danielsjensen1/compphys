import matplotlib.pyplot as plt
from numpy import array, linspace
from chap4.h2plus import H2Plus, Atom


#  Place the atoms one bohr radius apart
R1 = array([-.5e0, 0e0, 0e0])
R2 = -R1
Z = 1e0
exponents = array([13.00773, 1.962079, 0.444529, 0.1219492])
left = Atom(Z, exponents, R1)
right = Atom(Z, exponents, R2)
h2plus = H2Plus((left, right))
h2plus.variational()
def print_energies():
    for n in range(1, h2plus.bsize + 1):
        approx = h2plus.energy(n)
        print('E{0} = {1:+.7e}'.format(n, approx))
def plot_groundstate():
    x = linspace(-3, 3, 100)
    n = 1 #  Ground state
    approx = h2plus.eigfunc(n, x)
    plt.plot(x, approx, 'x-', label='STO-4G')
    plt.xlabel(r'$x$ (where $y=z=0$)')
    plt.ylabel(r'Ground state wave function')
    plt.legend()
    plt.title('Ground-state Electronic Eigenfunction of $H_{2}^{+}$ Molecule')
    plt.show()
#    plt.savefig('prob01Aplot.pdf')

if __name__ == '__main__':
    print_energies()
    plot_groundstate()
