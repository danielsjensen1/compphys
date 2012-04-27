import matplotlib.pyplot as plt
from numpy import append, arange, array, linspace, pi, piecewise, sqrt


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
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    ax1.plot(x, y)
    ax1.set_title('Kronig-Penny Potential')
    ax1.set_xticks(linspace(a/2e0, a/2e0+(barriers-1)*a, barriers))
    ax1.set_xticklabels(['${0}a$'.format(i) for i in range(barriers)])
    ax1.set_yticks([V0])
    ax1.set_yticklabels(['$V_0$'])
    ax1.set_ylim((0, V0+V0/2e0))
    plt.show()

def free_electron(a=2e0, Vmax=30e0, numpts=1000):
    kmax = 2 * pi / a
    E = linspace(0, Vmax, numpts)

    ax1 = plt.subplot2grid((1, 4), (0, 0), colspan=3)
    kplus = sqrt(2 * E)
    kminus = -kplus
    ax1.plot(kminus, E, 'r-', kplus, E, 'r-')
    zones = arange(kmax/2e0, kplus[-1], kmax)
    zones = append(-zones, zones)
    [plt.axvline(x=val, color='grey') for val in zones]
    ax1.set_xticks(zones)
    ax1.set_xticklabels(['${0:.1f}K$'.format(zone/kmax) for zone in zones])
    ax1.set_title('Extended-zone')
    
    ax2 = plt.subplot2grid((1, 4), (0, 3))
    kplus = (sqrt(2 * E) - kmax / 2e0) % kmax - kmax/2e0
    kminus = -kplus
    ax2.plot(kplus, E, 'r.', kminus, E, 'r.')
    ax2.set_xticks([-kmax/2e0, kmax/2e0])
    ax2.set_xticklabels(['$-1.5K$', '$1.5K$'])
    ax2.set_yticks([])
    ax2.set_title('Reduced Brillouin Zone')
    plt.show()

def muffin_tin(a1=array([1e0,0e0]), a2=array([0e0,1e0]), grid=(3,3)):
    pass

if __name__ == '__main__':
    plot_potential()
#    free_electron()