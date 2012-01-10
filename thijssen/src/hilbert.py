from numpy import dot, fromfunction
from scipy.linalg import eig


def hilbert(i, j):
    return ((i+1) + (j+1) - 1)**(-1)

N = 5
a = fromfunction(hilbert, (N, N), dtype=float)
print a
eigvals, eigvecs = eig(a)
print eigvals
print eigvecs
print dot(a, eigvecs[:,N-1]) / eigvecs[:,N-1]