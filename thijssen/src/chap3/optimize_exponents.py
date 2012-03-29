from numpy import array
from scipy.optimize import fmin
from scipy.optimize import fmin_l_bfgs_b
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
    guess = array([15e0, 10e0, 5e0, 1e0])
    guess = array([ 12.97135352,   1.94742281,   0.4415395 ,   0.12142323])
#    guess = array([12e0, 2e0, 1e0, .5e0])
    optimal_exponents = fmin_l_bfgs_b(opt_func, guess, approx_grad=True, 
                                      bounds=((1e-6, None),)*len(guess),
                                      epsilon=1e-11)
    print optimal_exponents
#    print energy
#    print func_calls