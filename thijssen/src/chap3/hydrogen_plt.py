import matplotlib.pyplot as plt
from numpy import array, linspace
from chap3.hydrogen import Hydrogen


exponents = array([13.00773, 1.962079, 0.444529, 0.1219492])
hydrogen = Hydrogen(exponents)
hydrogen.variational()
def print_energies():
    for n in range(1, len(exponents)+1):
        approx = hydrogen.approx_energy(n)
        exact = hydrogen.exact_energy(n)
        print('E{0} = {1:+.7e}, {2:.7e}'.format(n, approx, exact))
def plot_groundstate():
    r = linspace(0, 3, 100)
    n = 1 #  Ground state
    approx = hydrogen.approx_eigfunc(n, r)
    exact = hydrogen.exact_eigfunc(n, r)
    if approx[1] * exact[1] < 0:
        approx *= -1e0
    plt.plot(r, approx, 'r-', label='STO-4G')
    plt.plot(r, exact, 'b--', label='STO')
    plt.xlabel(r'$r$')
    plt.ylabel(r'$\phi_{1s}$')
    plt.legend()
    plt.title('Ground-state Electronic Eigenfunction of Hydrogen Atom')
#    plt.show()
    plt.savefig('prob03plot.pdf')

if __name__ == '__main__':
    print_energies()
    plot_groundstate()