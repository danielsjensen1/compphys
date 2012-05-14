from numpy import cos
from scipy.integrate import odeint


def duffing(t, x0, m=1e0, a=.25e0, b=0.5e0, F0=2.0e0, omega=2.4e0, gam=0.1e0):
    """
    Solution of the Duffing oscillator.
    
    Parameters
    ----------
    t : 1D array
        Time array containing the times when the solution should be returned.
    x0 : iterable
        Initial conditions such that x0[0] = initial position and
        x0[1] = initial velocity.
    """
    def func(y, t):
        x  = y[0]
        xp = y[1]
        xpp = (-gam * xp + 2 * a * x - 4 * b * x**3 + F0 * cos(omega * t)) / m
        return [xp, xpp]
    def gradient(y, t):
        x  = y[0]
        xp = y[1]
        return [[0e0, 1e0], [1/m*(2*a-12*b*x**2), -gam/m]]
    x = odeint(func, x0, t, Dfun=gradient, printmessg=True)
    return x