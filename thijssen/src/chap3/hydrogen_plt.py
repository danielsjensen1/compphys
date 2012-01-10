import matplotlib.pyplot as plt
from numpy import array, linspace
from chap3.hydrogen import Hydrogen


exponents = array([13.00773, 1.962079, 0.444529, 0.1219492])
hydrogen = Hydrogen(exponents)
hydrogen.variational()
r = linspace(0, 3, 100)
plt.plot(r, hydrogen.approx_eigfunc(0, r), 'r-', label='STO-4G')
plt.plot(r, hydrogen.exact_eigfunc(0, r), 'b--', label='STO')
plt.legend()
plt.show()