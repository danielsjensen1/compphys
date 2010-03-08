from numpy import array
from scipy.optimize import fmin
from chap3.hydrogen import Hydrogen


def opt_func(exponents):
    """Find the optimal exponents for the hydrogen wave function."""
    hydrogen = Hydrogen(exponents)
    hydrogen.variational()
    return hydrogen.eigvals[0]

if __name__ == '__main__':
    """You have to start fairly close to the exact answer to minimize this with
    the fmin routine or else it will start using negative or repeated values,
    which eigh can't handle."""
    guess = array([12e0, 2e0, 1e0, .5e0])
    print opt_func(guess)
    optimal_exponents = fmin(opt_func, guess)
    print optimal_exponents
#    print energy
#    print func_calls