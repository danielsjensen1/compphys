from numpy import linspace
from chap3.square_well import SquareWell
import matplotlib.pyplot as plt


sw = SquareWell()
#print sw.overlap(1, 2)
sw.variational()
x = linspace(-1, 1, 100)

for i in range(1, 3):#, sw.wfs):
    line1, = plt.plot(x, sw.exact_eigfunc(i+1, x), ls='-')
    clr = line1.get_color()
    line2, = plt.plot(x, sw.approx_eigfunc(i, x), color=clr, marker='.',
                      ls='None')
plt.title('Infinite Square Well')
plt.show()
explanation = """
The 'exact' solution is only exact up to an overall phase, hence the
difference in sign apparent in several of the eigenfunctions.  Notice that the
highest-order wave functions are the worst approximations while the lowest-order
wave functions are the best approximations.  Notice that there was no need to 
normalize our basis functions since the overlap matrix took care of that for us.

Also explain that the approximate basis functions have to satisfy the boundary
conditions.  I should probably plot the probability density instead of the wave
function to make comparisons easier.
"""
print explanation