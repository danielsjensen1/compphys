from numpy import cos, linspace
from scipy.integrate import odeint


m = 1e0
a = 0.25
b = 0.5
gam = 0.1
F0 = 1
omega = 2e0
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