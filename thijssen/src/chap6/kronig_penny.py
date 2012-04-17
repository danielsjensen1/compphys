from numpy import linspace, piecewise
import matplotlib.pyplot as plt

def plot_potential(a=2e0, Delta=1e0, V0=1.5e0, barriers=4, numpts=1000):
    """
    Plot the Kronig-Penny potential
    
    Parameters
    ----------
    a : real
        The direct-space lattice vector
    Delta : real
        The width of the potential.  This width must be less than `a`.
    V0 : real
        The height of the barrier
    """
    x = linspace(0e0, barriers*a, num=numpts)
    y = piecewise(x, [(x - Delta / 2e0) % a < Delta], [V0])
    plt.plot(x, y)
    plt.title('Kronig-Penny Potential')
    plt.show()

if __name__ == '__main__':
    plot_potential()