from numpy import indices, zeros


class DLA(object):
    def __init__(self, size=(175,175), occupied=None):
        #  Check size for odd integers
        self.size = size
        #  Create grid and make the center point occupied
        self.grid = zeros(size, dtype='bool')
        if occupied == None:
            middle = tuple(i / 2 for i in size)
            self.grid[middle] = True
        else:
            self.grid[occupied] = True
        ind = indices(self.grid.shape)
        self.grid[tuple(ind[:,:,-1])]
    
    def walker(self, pos):
        #  Add support for sticking coefficient
        pass
        