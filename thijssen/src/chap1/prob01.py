import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from numpy import arange, array, linspace, log, pi, poly1d, polyfit, where
from chap1.duffing import duffing


def partA(numpts=2000):
    t = linspace(0, 100, numpts)
    positions = (0.5e0, 0.5001e0)
    v0 = 0e0
    solutions = [duffing(t, [x0, v0], m=1e0, a=.5e0, b=0.25e0, F0=2.0e0, omega=2.4e0, gam=0.1e0) for x0 in positions]
    for pos, sol in zip(positions, solutions):
        plt.plot(t, sol[:,0], label='$x_0={0}$'.format(pos))
    plt.legend()
    plt.show()

def partB(numpts=2000):
    kwargs = {'F0':2e0, 'omega':2.4e0, 'gam':0.1e0, 'm':1e0, 'a':.5e0, 'b':.25e0}
    T = 2 * pi / kwargs['omega']
    t = arange(0, numpts*T, T)
#    t = linspace(0, 10*T, numpts)
    sol = duffing(t, x0=[0.5e0, 1e0], **kwargs)
    x, v = (sol[:, 0], sol[:, 1])
    plt.plot(x, kwargs['m']*v, color='black', marker=',', linestyle='...')
    plt.xlabel('$x$')
    plt.ylabel('$p$')
    plt.title('Strange Attractor for Duffing Oscillator')
    plt.show()

def partD(numpts=6400):
    kwargs = {'F0':2e0, 'omega':2.4e0, 'gam':0.1e0, 'm':1e0, 'a':.5e0, 'b':.25e0}
    T = 2 * pi / kwargs['omega']
    t = arange(0, numpts*T, T)
    sol = duffing(t, x0=[0.5e0, 1e0], **kwargs)
    x, p = (sol[:, 0], kwargs['m'] * sol[:, 1])
    xmin, xmax = (min(x), max(x))
    pmin, pmax = (min(p), max(p))
    Lx = xmax - xmin
    Lp = pmax - pmin
    L = max(Lx, Lp)
    def occupied(divisions):
        """
        Create array showing what squares are occupied.
        """
        N = float(divisions)
        inside = [[any((x >= xmin + i * L / N) * (x < xmin + (i + 1) * L / N) *
                       (p >= pmin + j * L / N) * (p < pmin + (j + 1) * L / N))
                   for j in range(divisions)]
                  for i in range(divisions)]
        return array(inside, dtype='bool')
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    attractor, = ax.plot(x, p, color='black', marker=',', linestyle='...')
    divisions = 2**4
    occ = occupied(divisions)
#    print(occ)
    print('sum(occ)={0}'.format(sum(sum(occ))))
    l = L / float(divisions)
    indices = zip(*where(occ))
    for (i, j) in indices:
        x_i = xmin + i * l
        p_i = pmin + j * l
        ax.add_patch(Rectangle((x_i, p_i), l, l, alpha=.1))
    plt.show()
    def count(divisions):
        """
        Find the number of occupied squares.
        
        This is a simpler and memory friendly version of `occupied` that only
        computes the total number of occupied squares.
        """
        N = float(divisions)
        occupied = sum(any((x >= xmin + i * L / N) * (x < xmin + (i + 1) * L / N) *
                           (p >= pmin + j * L / N) * (p < pmin + (j + 1) * L / N))
                       for i in range(divisions) for j in range(divisions))
        return occupied
    divisions = range(2, 8)
    logb = log(L / array(divisions, dtype='float'))
    num_occupied = map(count, divisions)
    print(num_occupied)
    lognum = log(num_occupied)
    z = polyfit(logb, lognum, 1)
    print(z)
    plt.plot(logb, lognum, 'r.')
    plt.plot(logb, poly1d(z)(logb), 'b-', label=r'$\log\left(N\right)={0:.3g}-{1:.3g}\cdot\log\left(b\right)$'.format(z[1], -z[0]))
    plt.legend()
    plt.show()


if __name__ == '__main__':
#    partA()
#    partB()
    partD()