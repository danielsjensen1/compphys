from itertools import product as iterproduct

import matplotlib.pyplot as plt
from numpy import any, array, asarray, min, product, zeros
from numpy.random import randint


class DLA(object):
    def __init__(self, size=(175, 175), occupied=None):
        self.size = size
        Ndim = len(size)
        #  Create grid and make the center point occupied
        self.grid = zeros(size, dtype=bool)
        if occupied is None:
            middle = tuple(i / 2 for i in size)
            self.grid[middle] = True
        else:
            self.grid[occupied] = True
        # Move to corners of hypercube
        self.directions = array(tuple(iterproduct((1, -1), repeat=Ndim)))
        # Move in only one direction at a time
#         self.directions = self._directions(Ndim)
        print(self.directions)
        self.Ndim = Ndim
    
    def _directions(self, Ndim):
        directions = zeros((2 * Ndim, Ndim), dtype=int)
        for i in xrange(Ndim):
            directions[i, i] = 1
            directions[-(i+1), i] = -1
        return directions
    
    def walker(self, initpos):
        grid, directions = self.grid, self.directions
        sizearr = array(self.size)
        Ndir = len(directions)
        #  TODO: Add support for sticking coefficient and history data.
        pos = asarray(initpos)
        stuck = False
        while(not stuck):
            direction = directions[randint(Ndir)]
            pos += direction
            pos = pos % sizearr
            neighbors = pos + directions
            if any([grid[tuple(neighbor % sizearr)] for neighbor in neighbors]):
                stuck = True
        grid[tuple(pos)] = True
        return pos
    
    def cluster(self, initpos, Nsites, Ninterval=200):
        Npts = product(self.size)
        logRg = []
        logN = []
        risum = 0e0
        risqrdsum = 0e0
        for i in xrange(min([Npts, Nsites - 1])):
            pos = self.walker(initpos)
            risum += pos
            risqrdsum += pos.dot(pos)
            if i % Ninterval:
                logRg = log()
        return self.grid
    
    def plt_lattice(self):
        fig, ax0 = plt.subplots(1, 1)
        ax0.matshow(self.grid)
        plt.show()