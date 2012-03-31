from numpy import array
from scipy.optimize import fmin
from scipy.optimize import fmin_l_bfgs_b
from scipy.optimize import fmin_slsqp
from chap3.hydrogen import Hydrogen


def opt_func(exponents):
    """Find the optimal exponents for the hydrogen wave function."""
#    print('exponents=', exponents)
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
#    guess = array([15e0, 10e0, 5e0, 1e0])
#    guess = array([12e0, 2e0, 1e0, .5e0])
    optimal_exponents = fmin_l_bfgs_b(opt_func, guess, approx_grad=True, 
                                      bounds=((1e-6, None),)*len(guess),
                                      epsilon=1e-11)
    print optimal_exponents
    guess = array([2**n for n in range(len(guess))[::-1]])
    print('guess=',guess)
    def ieqcons(x):
        tol = 1e-1
        return array([x[0]-x[1]-tol, x[1]-x[2]-tol, x[2]-x[3]-tol, x[3]-tol])
    def prime_ieqcons(x):
        return array([[1e0, -1e0,  0e0,  0e0],
                      [0e0,  1e0, -1e0,  0e0],
                      [0e0,  0e0,  1e0, -1e0],
                      [0e0,  0e0,  0e0,  1e0]])
    optimal_exponents = fmin_slsqp(opt_func, guess, f_ieqcons=ieqcons, 
                                   fprime_ieqcons=prime_ieqcons,
                                   acc=1e-11)
    print optimal_exponents
    
    print('"Exact"')
    guess = array([13.00773, 1.962079, 0.444529, 0.1219492])
    print opt_func(guess)