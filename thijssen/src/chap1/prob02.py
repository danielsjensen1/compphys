from chap1.dla import DLA


def partA(Nsites=9000, size=(150, 150)):
    dla = DLA(size=size, occupied=None)
    dla.cluster(initpos=(0, 0), Nsites=Nsites)
    dla.plt_lattice()


if __name__ == '__main__':
    partA()