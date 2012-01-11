import matplotlib.pyplot as plt
from numpy import cos, linspace, pi
from scipy.integrate import odeint


class Duffing(object):
    def __init__(self, m=1e0, a=.25e0, b=0.5e0):
        pass

gam = 0.1
F0 = 1
omega = 2e0
T = 2 * pi / omega
x0 = [0e0, 1e0]
#  Part a
def func(y, t):
    x  = y[0]
    xp = y[1]
    xpp = (-gam * xp + 2 * a * x - 4 * b * x**3 + F0*cos(omega * t)) / m
    return [xp, xpp]

numpts = 100
t = linspace(0, 10, numpts)
x = odeint(func, x0, t)
print x.shape
plt.plot(t, x[:,0], label='position')
plt.show()