from numpy import linspace
import matplotlib.pyplot as plt
from chap3.prob03squarewell import SquareWell


def make_plots(N):
    sw = SquareWell(N)
    sw.variational()
    approx_eigvals = sw.eigvals
    x = linspace(-1, 1, 100)
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(1, 1, 1)
    for i in range(0, 5):#, sw.wfs):
        exact = sw.exact_eigfunc(i+1, x)
        approx = sw.approx_eigfunc(i, x)
        if exact[1] * approx[1] < 0:
            approx = approx * -1e0
        line1, = ax.plot(x, exact, ls='-')
        clr = line1.get_color()    
        line2, = ax.plot(x, approx, color=clr, marker='.',
                          ls='None', label='{0:.5g} hartrees'.format(approx_eigvals[i]))
    ax.set_xlabel('x')
    ax.set_ylim(-1, 1)
    plt.title('Infinite Square Well')
    plt.legend()
    filename = 'image{0}.pdf'.format(N)
    plt.savefig(filename)

if __name__ == '__main__':
    Ns = (5, 8, 12, 16)
    for N in Ns:
        make_plots(N)