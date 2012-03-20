import itertools
import matplotlib.pyplot as plt
from numpy import arange, argmin, array, exp, linspace
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
    E = array([h2plus_calc(d, exponents).energy_tot(n) for d in darray])
#    print min(E)
#    print darray[argmin(E)]
    return E

def pes_plots(darray):
    markers = itertools.cycle(('+', '.', '-', '|', '*', ','))
    for name, exponents in sorted(H_exponents.items()):
        energies = pes(darray, exponents)
        plt.plot(darray, energies, label=name, marker=markers.next(), linewidth=.5)
        print('Approximation = {0}, Energy = {1:.5e} at d = {2:.5g}'.format(name, min(energies), darray[argmin(energies)]))
    plt.xlabel(r'Interatomic distance $d$')
    plt.xlim((1e0, 3e0))
    plt.ylabel(r'Ground state energy (Hartrees)')
    plt.legend()
    plt.title('Potential Energy Curve(s) for $H_{2}^{+}$')
    plt.show()
#    plt.savefig('prob01Aplot.pdf')

def plot_basis(rarray, norm=True):
    i = 0
    fig = plt.figure(figsize=(7, 10))
    if norm:
        title = "Basis Functions (norm=1)"
    else:
        title = "Basis Functions (max=1)"
    fig.suptitle(title, fontsize=18)
    for name, exponents in sorted(H_exponents.items()):
        i += 1
        markers = itertools.cycle(('+', '.', '-', '|', '*', ','))
        for alpha in exponents:
            if norm:
                N = (2 * alpha / 4e0) ** (0.75e0) #  Normalization constant
            else:
                N = 1e0
            phi = N * exp(-alpha * rarray ** 2)
            plt.subplot(3, 2, i)
            plt.plot(rarray, phi, dashes=(2,2),label='{0:.5g}'.format(alpha), marker=markers.next(), linewidth=.5)
            leg = plt.legend(numpoints=1)
            for t in leg.get_texts():
                t.set_fontsize('small')    # the legend text fontsize
            for l in leg.get_lines():
                l.set_linewidth(0.5)  # the legend line width
            plt.xlabel(r'$r$')
            plt.ylabel(r'Basis function')
            plt.title(name)
    fig.subplots_adjust(hspace=0.4, wspace=0.3)
    plt.show()
#    if norm:
#        plt.savefig('basis_normed.pdf')
#    else:
#        plt.savefig('basis_maxed.pdf')

if __name__ == '__main__':
    darray = arange(1e0, 3.1e0, 0.01)
    pes_plots(darray)

#    exponents = H_exponents['STO-4G']
#    pes(darray, exponents)

    rarray = linspace(0, 2e0, 50)
    plot_basis(rarray)
    plot_basis(rarray, norm=False)
